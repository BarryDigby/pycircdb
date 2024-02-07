from pyveplot import Hiveplot, Axis, Node
import random

# Create a Hiveplot object
h = Hiveplot()

# Create 3 axes
h.axes = [Axis(start=20, angle=0,
                stroke=random.choice(c), stroke_width=1.1),
            Axis(start=20, angle=120,
                stroke=random.choice(c), stroke_width=1.1),
            Axis(start=20, angle=120 + 120,
                stroke=random.choice(c), stroke_width=1.1)
            ]

# Create 3 nodes on each axis
for i in range(3):
    circRNAs_node = Node(i, random.uniform(0.1, 0.9), stroke='blue', fill='blue')
    miRNA_node = Node(i, random.uniform(0.1, 0.9), stroke='red', fill='red')
    mRNAs_node = Node(i, random.uniform(0.1, 0.9), stroke='green', fill='green')

    # Add the nodes to the axes
    circRNAs_axis.add_node(circRNAs_node)
    miRNA_axis.add_node(miRNA_node)
    mRNAs_axis.add_node(mRNAs_node)

    # Add edges between the nodes
    if i > 0:
        h.connect(circRNAs_axis.nodes[i-1], miRNA_axis.nodes[i])
        h.connect(miRNA_axis.nodes[i-1], mRNAs_axis.nodes[i])
        h.connect(mRNAs_axis.nodes[i-1], circRNAs_axis.nodes[i])

# Save the Hiveplot to a file
h.save('foo.svg')