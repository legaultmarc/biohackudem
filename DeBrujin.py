class DeBruijnGraph(object):
    """ A de bruijn multigraph built from a collection of strings.
    """

    reverse_rule = {'C':['C'], 'G':['G', 'A'], 'T':['C', 'T'], 'A':'A'}

    @staticmethod
    def toKmer(strg, k):
        """return all the kmer of len k found in a String
        """
        for i in xrange(0, len(strg)-(k-1)):
            yield strg[i:i+k]

    @staticmethod
    def node_from_edge(edge_seq):
        k = len(edge_seq)
        return edge_seq[0:k-1], edge_seq[1:k]

    class Edge(object):
        """An edge that link two node
        """
        def __init__(self, prev, next, chev):
            self.prev = prev
            self.next = next
            self.chev = chev
            self.visited = 0;

        def __repr__(self):
            return "%s--(%s)--%s"%(self.prev, self.chev, self.next)

        def __hash__(self):
            return hash(self.prev)+hash(self.chev)+hash(self.next)
        def __eq__(self, other):
            return isinstance(other, self.__class__) and\
                    self.chev==other.chev and\
                    self.prev==other.prev and\
                    self.next==other.next

    class Node(object):
        """Add a node in a De bruijn graph"""
        def __init__(self, kmerseq="$", color=1):
            self.kmerseq = kmerseq
            self.edge_list = set()
            self.color = color
            # keeping this will make fetching for adjacent node easier
            self.adjacent_nodes = set()

        def add_edge(self, next_node, chev):
            self.edge_list.add(DeBruijnGraph.Edge(self, next_node, chev))
            self.adjacent_nodes.add(next_node)

        def __hash__(self):
            return hash(self.kmerseq)

        def __repr__(self):
            return "%s(%i)"%(self.kmerseq, self.color)

        def __eq__(self, other):
            # Comparision of node
            # We just compare the sequence
            return isinstance(other, self.__class__) and self.kmerseq==other.kmerseq

        def is_root(self):
            return self.kmerseq=="$"

        def is_adjacent(self, search_node):
            return next((node for node in self.adjacent_nodes if node==search_node), None)

        def has_adjacent(self):
            return len(self.adjacent_nodes) > 0

        def get_possible_chev(self, othernode):
            return self.kmerseq+othernode.kmer[-1]

    def __init__(self, readlist, k, colorcode=1):
        # this graph contains the list of kmer and the corresponding note
        self.graph = {}
        self.k = k

        for read in readlist:
            for edge in self.toKmer(read, k):
                # add a new node or return the node at key in the dict
                previous_node, next_node = self.node_from_edge(edge)
                self.graph[previous_node] = self.graph.get(previous_node, DeBruijnGraph.Node(previous_node, color=colorcode))
                self.graph[next_node] = self.graph.get(next_node, DeBruijnGraph.Node(next_node, color=colorcode))
                self.graph[previous_node].add_edge(self.graph[next_node], edge)

    def bfs(self):
        visited, queue = [], self.graph.values()
        while queue:
            vertex = queue.pop(0)
            if vertex not in visited:
                visited.append(vertex)
                queue.extend(vertex.adjacent_nodes - set(visited))

        return visited

    def find_polymorphisme(self):
        pass

if __name__ == '__main__':

    strlist = ['ATGGCGTGCAATG']

    graph = DeBruijnGraph(strlist, 3)
    for key, value in graph.graph.items():
        print "%s : %s"%(value, value.edge_list)
    print "\n\n"
    print graph.bfs()
