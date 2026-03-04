import pandas as pd
import networkx as nx
import streamlit as st
from streamlit_d3graph import d3graph


@st.cache_data
def init_graph():
    # Initialize
    d3 = d3graph()
    # Load my example
    df = pd.read_csv("results/filtered_mirnas.csv")

    # Create Graph
    G = nx.from_pandas_edgelist(df, "hg38", "miRNA", create_using=nx.DiGraph())
    print(G.nodes())   
    print(len(G.nodes()))

    # Adj mat 
    adjmat = pd.DataFrame(nx.adjacency_matrix(G).toarray())

    # Calculate the degrees for each node
    degrees = [deg for node, deg in G.degree()]

    # Create a DataFrame with the degrees and the node labels
    df_degrees = pd.DataFrame({
        'degree': degrees,
        'label': G.nodes()
    })

    print(df_degrees)

    label = df_degrees["label"].values
    node_size = df_degrees["degree"].values

    return d3, adjmat, df_degrees, label, node_size


@st.cache_data
def graph_one(_d3, adjmat, df, label, node_size):
    d3.graph(adjmat)
    d3.set_node_properties(color=df["label"].values)
    return d3


@st.cache_data
def graph_two(_d3, adjmat, df, label, node_size):
    d3.graph(adjmat)
    d3.set_node_properties(label=label, color=label, cmap="Set1")
    return d3


@st.cache_data
def graph_three(_d3, adjmat, df, label, node_size):
    d3.set_node_properties(color=label, size=node_size)
    return d3


@st.cache_data
def graph_four(_d3, adjmat, df, label, node_size):
    d3.set_edge_properties(edge_distance=100)
    d3.set_node_properties(color=node_size, size=node_size)
    return d3


@st.cache_data
def graph_five(adjmat, df, label, node_size):
    d3 = d3graph(charge=1000)
    d3.graph(adjmat)
    d3.set_node_properties(color=node_size, size=node_size)
    return d3


@st.cache_data
def graph_six(adjmat, df, label, node_size):
    d3 = d3graph(collision=1, charge=250)
    d3.graph(adjmat)
    d3.set_node_properties(
        color=label, size=node_size, edge_size=node_size, cmap="Set1"
    )
    return d3


@st.cache_data
def graph_seven(adjmat, df, label, node_size):
    d3 = d3graph(collision=1, charge=250)
    d3.graph(adjmat)
    d3.set_node_properties(
        color=label,
        size=node_size,
        edge_size=node_size,
        edge_color="#00FFFF",
        cmap="Set1",
    )
    return d3


d3, adjmat, df, label, node_size = init_graph()

#d3 = graph_one(d3, adjmat, df, label, node_size)
#d3.show()

d3 = graph_two(d3, adjmat, df, label, node_size)
d3.show()

#d3 = graph_three(d3, adjmat, df, label, node_size)
#d3.show()
#
#d3 = graph_four(d3, adjmat, df, label, node_size)
#d3.show()
#
#d3 = graph_five(adjmat, df, label, node_size)
#d3.show()
#
#d3 = graph_six(adjmat, df, label, node_size)
#d3.show()
#
#d3 = graph_seven(adjmat, df, label, node_size)
#d3.show()