# -*- coding: utf-8 -*-
"""
Created by Gaurang Mahajan on June 21 2021.
"""

import streamlit as st
import numpy as np
import pandas as pd
from scipy import stats
from operator import itemgetter
import networkx as nx


def targets(tf,g):
    return [s for s in g.successors(tf)]

def find_best_tf(tfs,bkg,sig,trn):
    tflist=[]
    for tf in tfs:
        a=len(np.intersect1d(targets(tf,trn),sig))
        b=len(sig)-a
        c=len(np.intersect1d(targets(tf,trn),bkg))-a
        d=len(bkg)-a-b-c
        _, p = stats.fisher_exact([[a,b],[c,d]],alternative='greater')
        tflist.append((tf,p))
    return sorted(tflist,key=itemgetter(1))[0][0]

def combined_pval(tfsol,allgenes,siggenes,trn):
    l = [set(targets(tf,trn)) for tf in tfsol]
    targetset = list(frozenset().union(*l))
    a=len(np.intersect1d(targetset,siggenes))
    b=len(siggenes)-a
    c=len(np.intersect1d(targetset,allgenes))-a
    d=len(allgenes)-a-b-c
    _, p = stats.fisher_exact([[a,b],[c,d]],alternative='greater')   
    return np.log10(p)
    

with st.sidebar:
    
    whichdata = st.radio("Run tool for:",("User input","Example"),1)

    if whichdata=='Example':
        st.markdown('Run on example (*E. coli*, microarray expression data for [pH 8.7 vs. pH 7](https://doi.org/10.1128/jb.187.1.304-319.2005) comparison, [RegulonDB 8.6](https://doi.org/10.1093/nar/gky1077))')
        network_file = './data/RegulonDB_pairs.txt'
        geneexp_file = './data/ph8_vs_ph7_degenes_fdr1e-3_fc2.txt'
    
    if whichdata=='User input':
        st.header("File uploads")
        gex_file = st.file_uploader("1. Upload tab-delimited diff. gene expression data (Format: 'Gene -> DGE status (0/1)', one per line)",type=['txt'])
        net_file = st.file_uploader("2. Upload tab-delimited regulatory network file (Format: 'TF -> gene', one per line)",type=['txt'])
        if gex_file is not None: geneexp_file = gex_file
        if net_file is not None: network_file = net_file
        
    st.header("Input parameters")
    with st.form(key="grid_reset"):
        n_cutoff = st.slider("Cutoff on TF network degree", 2, 20)
        p_cutoff = st.number_input("-Log p-val cutoff\n on indiv. TF enrichment", 2, 6)
        #chart_type = st.selectbox("Choose chart type", plot_types)
        submit_status = st.form_submit_button(label="Score TFs on input dataset")
    
    with st.expander("About this app"):
        st.markdown("Computes differentially active regulators (trans factors/TF), given a user-specified 'large' list of gene expression changes between two conditions\
        and a precompiled TF-gene regulatory network (implements the approach proposed in \
        [Mahajan & Mande, 2015](https://doi.org/10.1371/journal.pone.0142147))")

with st.container():
    st.header("**Relating genome-wide expression changes to altered regulatory protein activity (Loaded DGE data -> TF scores)**")
   
e1=0
e2=0

if submit_status:
    ####################################
    #### Reading network file and creating a Networkx directed network object:

    g = nx.DiGraph()
    try:
        trn = pd.read_csv(network_file,header=None,delimiter='\t')
        for s,t in zip(trn[0],trn[1]): g.add_edge(s,t)
    except IOError as e1:
        with st.container():
            st.write("I/O error({0}): {1}".format(e1.errno, e1))

    regs = [r for r in g.nodes() if g.out_degree(r) >= n_cutoff]
    with st.container(): 
        st.write('Pre-compiled regulatory network comprises',g.number_of_nodes(),'genes,',\
          g.number_of_edges(),'directed links, and',len(regs),'TFs with out-degree >',n_cutoff)
    ####################################

    ####################################
    #### Reading in list of all expressed genes and significantly DE genes:

    siggenes=[]
    allgenes=[]
    try:
        gex = pd.read_csv(geneexp_file,header=None,delimiter='\t')
        for gid,e in zip(gex[0],gex[1]):
            allgenes.append(gid)
            if float(e)==1: siggenes.append(gid)
    except IOError as e2:
        with st.container():
            st.write("I/O error({0}): {1}".format(e2.errno, e2))
        
    siggenes=list(set(siggenes))
    allgenes=list(set(allgenes))
    tflist = [tf for tf in regs if len(np.intersect1d(siggenes,targets(tf,g)))>0]
    with st.container():
        st.write(len(siggenes),'genes out of a total of', len(allgenes),'are labeled as significantly diff. expressed (DEG).')
        st.write('Target sets of',len(tflist),'TFs have non-zero overlap with the DEG set.')
    ####################################

    if e1!=0 or e2!=0 or len(tflist)==0:
        with st.container(): st.write("\nError: Can't run test.")
        
    if ((e1==0) and (e2==0) and (len(tflist)>0)):

        ####################################
        #### FET enrichment p-values for individual TFs:

        sigtfs={}
        for tf in tflist[:]:
           a = len(np.intersect1d(siggenes,targets(tf,g)))
           b = len(siggenes)-a
           c = len(np.intersect1d(allgenes,targets(tf,g)))-a
           d = len(allgenes)-a-b-c
           o,p = stats.fisher_exact([[a,b],[c,d]],alternative='greater')
           sigtfs[tf] = np.min([1,p*len(tflist)])
        with st.container(): 
            st.write((len([tf for tf in sigtfs if sigtfs[tf] < 10**(-p_cutoff)])),\
            ' TFs are individually associated with DEGs at p(adj)<',10**(-p_cutoff),'significance level:')
            
        df1 = pd.DataFrame(sorted([[tf,len(targets(tf,g)),len(np.intersect1d(siggenes,targets(tf,g))),sigtfs[tf]] for tf in sigtfs \
        if sigtfs[tf] < 10**(-p_cutoff)],key=itemgetter(-1)),\
        columns=['TF','Total # of targets','# DEG targets','Indiv. p-val\n(FET; adjusted)'])
        with st.container(): df1
        ####################################

        with st.container(): 
            st.markdown("**TF subnetwork output**<sup>a</sup>:",unsafe_allow_html=True)
            
        #########################################
        #### Greedy approach for finding the optimal TF combination ('Method C' in cited paper):

        tfset=tflist[:]
        tf = find_best_tf(tfset,allgenes,siggenes,g)
        tfsol_old=[tf]
        tfset.remove(tf)
        tfset_rem=tfset[:]
        old_pval=combined_pval(tfsol_old,allgenes,siggenes,g)
        if len(tfset)>0:
           while len(tfset_rem[:])>0: 
                 tflist=[]
                 for tf in tfset_rem:
                     tfsol_temp = tfsol_old[:]
                     tfsol_temp.append(tf)
                     tflist.append((tf,combined_pval(tfsol_temp,allgenes,siggenes,g)))
                 (tf,p)=sorted(tflist,key=itemgetter(1))[0]
                 if p >= old_pval:
                    break
                 else:  
                    tfsol_new=tfsol_old[:]
                    tfsol_new.append(tf)
                    old_pval = p
                    tfsol_old = tfsol_new[:]
                    tfset_rem = [r for r in tfset_rem[:] if r!=tf]
        tfsol_opt = tfsol_old[:]
        pval_opt = combined_pval(tfsol_opt,allgenes,siggenes,g)
        
        df2 = pd.DataFrame(sorted([[tf,len(targets(tf,g)),len(np.intersect1d(siggenes,targets(tf,g))),sigtfs[tf]] for tf in tfsol_opt],\
        key=itemgetter(-1)),columns=['TF','Total # of targets','# DEG targets','Indiv. p-val\n(FET; adjusted)'])
        with st.container(): 
            st.write('Optimal aggregated p-value (log10) =',round(pval_opt,1))
            df2
        
        ###############################################    
        #### Null distribution generated by repeated randomization of input DGE data, and running the greedy approach on each instance:

        with st.container(): 
            st.markdown('**Distribution of min. p-values on randomized data (*-ve control* based on 20 randomized instances) <sup>b</sup>**:',unsafe_allow_html=True)           
        pval_dist = []
        tfscore = {tf:0 for tf in regs}
        Nmax = 20
        for n in range(Nmax):
        
            siggenes_rnd=list(np.random.choice(allgenes,len(siggenes)))
            tflist = [tf for tf in regs if len(np.intersect1d(siggenes_rnd,targets(tf,g)))>0]
        
            tfset=tflist[:]
            tf = find_best_tf(tfset,allgenes,siggenes_rnd,g)
            tfsol_old=[tf]
            tfset.remove(tf)
            tfset_rem=tfset[:]
            old_pval=combined_pval(tfsol_old,allgenes,siggenes_rnd,g)
            if len(tfset)>0:
                while len(tfset_rem[:])>0: 
                    tflist=[]
                    for tf in tfset_rem:
                        tfsol_temp = tfsol_old[:]
                        tfsol_temp.append(tf)
                        tflist.append((tf,combined_pval(tfsol_temp,allgenes,siggenes_rnd,g)))
                    (tf,p)=sorted(tflist,key=itemgetter(1))[0]
                    if p >= old_pval:
                        break
                    else:  
                        tfsol_new=tfsol_old[:]
                        tfsol_new.append(tf)
                        old_pval = p
                        tfsol_old = tfsol_new[:]
                        tfset_rem = [r for r in tfset_rem[:] if r!=tf]
            tfsol_approx_rnd = tfsol_old[:]
            pval_rnd = combined_pval(tfsol_approx_rnd,allgenes,siggenes_rnd,g)
            pval_dist.append(pval_rnd)
            for tf in tfsol_approx_rnd: tfscore[tf] +=1
         
        df3 = sorted([(tf,100*tfscore[tf]/float(Nmax)) for tf in tfsol_opt],key=itemgetter(-1),reverse=True)
        st.write('z-score of the original TF subnet output = ',round((pval_opt - np.mean(pval_dist))/np.std(pval_dist),1))
        
        ######################################
        #### Plotting the summary statistics:

        fig = plt.figure(figsize=(6,2))
        
        ax = fig.add_subplot(131)
        plt.axvline(pval_opt,color='gray',linestyle='--')
        ax.hist(pval_dist,bins=10,alpha=0.7)
        ax.set_xlabel('Log10 (p-val)')
        ax.set_ylabel('Frequency')
        ax = fig.add_subplot(133)
        ax.barh(range(1,1+len(tfsol_opt)),[s for t,s in df3],height=0.3,align='center')
        ax.set_yticks(range(1,1+len(tfsol_opt)))
        ax.set_yticklabels([t for t,s in df3])
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Occurrence % in\nrandom trials')
        ax.set_ylabel('TFs in best-scoring combination')
        ax.set_xlim([0,100])
        
        st.pyplot(fig)
        #####################################
        
    st.markdown("***")
    st.markdown("Notes:<br />1) <sup>a</sup>Greedy bottom-up search is used to estimate the TF combination, the aggregated target set of\
 which has minimum p-value (one-sided Fisher's exact test) for association with the list of DEGs.\
<br />2) <sup>b</sup>The obtained match (log p-value) is compared with a null distribution which is generated by repeated randomization of the original DEG\
 set. A large deviation from this distribution (large -ve z-score) indicates a concordance between the regulatory network\
 and the DEG set, and is suggestive of a potential set of hypotheses (regulatory subnetwork) to account for the differential expression profile.\
 <br />3) Raw p-values for individual regulators are adjusted (Bonferroni correction) to account for multiple testing.\
<br />4) More details about the method may be found [here](https://doi.org/10.1371/journal.pone.0142147).\
<br />5) The uploaded data is not stored by the app.",unsafe_allow_html=True)
    st.markdown("***")
    st.markdown("Created by Gaurang Mahajan. Any [feedback](https://gaurangm-web.github.io) on use of this tool or suggestions for improvement will be appreciated.",unsafe_allow_html=True)