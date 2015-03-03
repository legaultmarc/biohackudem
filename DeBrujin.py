class DeBruijnGraph(object):
    """ A de bruijn multigraph built from a collection of strings.
    """

    reverse_rule = {'C':['C'], 'G':['G', 'A'], 'T':['C', 'T'], 'A':'A'}

    @staticmethod
    def toKmer(strg, k, color=1):
        """return all the kmer of len k found in a String
        """
        for i in xrange(0, len(strg)-(k-1)):
            yield DeBruijnGraph.Node(strg[i:i+k], color=color)

    class Edge(object):
        """An edge that link two node
        """
        def __init__(self, prev, next, chev):
            self.prev = prev
            self.next = next
            self.chev = chev

    class Node(object):
        """Add a node in a De bruijn graph"""
        def __init__(self, kmerseq="$", color=1):
            self.kmerseq = kmerseq
            self.edge_list = []
            self.visited = 0;
            self.color = color
            # keeping this will make fetching for adjacent node easier
            self.adjacent_nodes = set()

        def add_edge(self, next_node, chev):
            self.edge_list.append(DeBruijnGraph.Edge(self, next_node, chev))
            self.adjacent_nodes.add(next_node)

        def __hash__(self):
            return hash(self.kmerseq)

        def __repr__(self):
            return "%s:%s(%i)"%(self.kmerseq, self.visited, self.color)

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

        def get_possible_chev(self):
            return self.kmerseq[-1]

    def __init__(self, readlist, k, colorcode=1):
        # this graph contains the list of kmer and the corresponding note
        self.graph = {}
        self.k = k

        for read in readlist:
            previous_node = None
            for node in self.toKmer(read, k, colorcode):
                # add a new node or return the node at key in the dict
                self.graph[node.kmerseq] = self.graph.get(node.kmerseq, node)
                self.graph[node.kmerseq].visited+=1

                # this is the first node
                if not previous_node:
                    previous_node = node

                # we have a new node, add edge between new node and previous node
                else :
                    self.graph[previous_node.kmerseq].add_edge(self.graph[node.kmerseq], self.graph[node.kmerseq].get_possible_chev())


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

    strlist = ['CGATTCTAAGT', 'CGATTGTAAGT']

    graph = DeBruijnGraph(strlist, 3)
    print graph.bfs()
