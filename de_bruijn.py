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

        self.value = seq[:k] if len(seq) > k else seq

        # Memorize this node.
        node_pointers[self.value] = self

        # Add adjacent.
        if len(seq) != k:
            self.edges.append((seq[k], Node.node_factory(seq[1:], k)))

    @classmethod
    def node_factory(cls, seq, k):
        if len(seq) == k:
            return cls(seq, k)
        elif seq[:k] in node_pointers:
            return node_pointers[seq[:k]]
        else:
            return cls(seq, k)

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
