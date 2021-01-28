# DiscreteVectorFieldCubicalComplex.py
# MIT LICENSE 2020 Ewerton R. Vieira

import pychomp

from DSGRN import DirectedAcyclicGraph


class DiscreteVectorFieldCubicalComplex:

    def complex(self):
        return self.cc

    def diagram(self):
        return self.digraph.adjacencies

    def closedrightfringe(self, i):
        return any(self.cc.rightfringe(k) for k in self.cc.star({i}))

    def blowupfringe(self, i):
        """
        determine if cell i (in dual complex) is in the fringe
        """
        return all(self.closedrightfringe(k) for k in self.oc.simplex(self.fc_to_cc(i)))

    def blowupinfinity(self, i):
        """
        determine if cell i (in cubical complex) is in the fringe,
        i. e., associate with the repeller at infinity
        """
        return self.blowupfringe(self.cc_to_fc(i))

    def cc_open(self):  # return cells that are in interior of cubical complex (not fringe cells/wrap cells)
        cc_open = []
        for j in range(self.cc.size()):
            if not self.closedrightfringe(j):
                cc_open.append(j)
        return cc_open

    def cc_closure(self):  # return cells that are fringe cells/wrap cells
        cc_closure = [a for a in self.cc]
        return list(set(cc_closure) - set(self.cc_open()))

    def top_cc_closure(self):  # return top cells that are fringe cells/wrap cells
        return list(set(self.cc_closure()).intersection(self.cc(self.D)))

    # return top cells that are in interior of cubical complex (not fringe cells/wrap cells)
    def top_cc_open(self):
        return list(set(self.cc_open()).intersection(self.cc(self.D)))

    # add s -> t in combinatorial complex, have to use cc.top_cell (true one) indexation
    def add_pair(self, s, t):
        self.digraph.add_edge(s, t)

    def remove_pair(self, s, t):
        self.digraph.remove_edge(s, t)

    def add_pairs(self, Add):  # add s : { } in combinatorial complex
        # using 0 to ? indexation for non fringe top cells (easy to write by hand)
        if list(Add.keys())[0] < self.cc.count()[0]:
            A = self.top_cc_open()
            A.sort()
            for _, a in enumerate(Add.keys()):
                for _, b in enumerate(Add.get(a)):
                    self.digraph.add_edge(A[a], A[b])

        else:  # using cc.top_cell (true one) indexation
            for _, a in enumerate(Add.keys()):
                for _, b in enumerate(Add.get(a)):
                    self.digraph.add_edge(a, b)

    def remove_pairs(self, Add):  # add s : { } in combinatorial complex
        # using 0 to ? indexation for non fringe top cells (easy to write by hand)
        if list(Add.keys())[0] < self.cc.count()[0]:
            A = self.top_cc_open()
            A.sort()
            for _, a in enumerate(Add.keys()):
                for _, b in enumerate(Add.get(a)):
                    self.digraph.remove_edge(A[a], A[b])

        else:  # using cc.top_cell (true one) indexation
            for _, a in enumerate(Add.keys()):
                for _, b in enumerate(Add.get(a)):
                    self.digraph.remove_edge(a, b)

    # return the adjacent top cell in direction d of a given top cell s

    def adjacent_top_cell_right(self, s, d):
        return self.cc.right(self.cc.right(s, d), d)

    # return the adjacent top cell in direction d of a given top cell s
    def adjacent_top_cell_left(self, s, d):
        return self.cc.left(self.cc.left(s, d), d)

    # put double edges for top cells in cc_closure
    def combinatorial_fringe(self):
        for s in self.top_cc_closure():
            for j in range(self.D):
                r_j = self.adjacent_top_cell_right(s, j)
                l_j = self.adjacent_top_cell_left(s, j)

                if r_j in self.top_cc_closure():
                    self.digraph.add_edge(s, r_j)
                    self.digraph.add_edge(r_j, s)
                else:
                    if self.att:
                        self.digraph.add_edge(s, r_j)
                    else:
                        self.digraph.add_edge(r_j, s)

                if l_j in self.top_cc_closure():
                    self.digraph.add_edge(s, l_j)
                    self.digraph.add_edge(l_j, s)
                else:
                    if self.att:
                        self.digraph.add_edge(s, l_j)
                    else:
                        self.digraph.add_edge(l_j, s)

    def auto_double_edges(self):  # add s : { } in combinatorial complex
        for s in self.top_cc_open():
            for j in range(self.D):
                r_j = self.adjacent_top_cell_right(s, j)
                l_j = self.adjacent_top_cell_left(s, j)

                if all([(s, r_j) not in self.digraph.edges(), (r_j, s) not in self.digraph.edges()]):
                    self.digraph.add_edge(s, r_j)
                    self.digraph.add_edge(r_j, s)

                if all([(s, l_j) not in self.digraph.edges(), (l_j, s) not in self.digraph.edges()]):
                    self.digraph.add_edge(s, l_j)
                    self.digraph.add_edge(l_j, s)

    def __init__(self, cc, att=True):
        self.cc = cc  # cubical complex
        self.D = cc.dimension()
        self.att = att  # true for attractor boundary
        self.digraph = DirectedAcyclicGraph()
        for cell in self.cc(self.D):  # combinatorial dynamical system
            self.digraph.add_vertex(cell)

        self.combinatorial_fringe()
