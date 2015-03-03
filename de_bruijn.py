#!/usr/bin/env python

import os

from Bio import SeqIO

DATA = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "data"
))


node_pointers = {}


class Node():
    def __init__(self, seq, k):
        # Tuples of right-character to node.
        self.edges = []

        # The node value is the kmer.
        self.value = seq[:k] if len(seq) > k else seq

        # Add adjacent.
        if len(seq) != k:
            self.edges.append((seq[k], Node(seq[1:], k)))

    def __repr__(self):
        return "<{}>".format(self.value)


def print_graph(g, depth=0):
    """Print a graph."""
    print " " * depth + repr(g)

    for lchar, node in g.edges:
        print_graph(node, depth + 1)
        
if __name__ == "__main__":
    k = 3
    g = None
    with open(os.path.join(DATA, "all.fastq.1")) as f:
        for record in SeqIO.parse(f, "fastq"):
            if g is None:
                g = Node("BANANA", k)
            else:
                pass
                #g.add_sequence(record.seq)

    print_graph(g)
