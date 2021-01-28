# DrawPosetGraph.py  # 2020-09-12
# MIT LICENSE 2020 Ewerton R. Vieira & Marcio Gameiro


from collections import Counter
import graphviz
import matplotlib.colors as mcolors
from colour import Color

class DrawPosetGraph():

    def PositionInList(self,j):#Given a list of coloring P and a vertex j return the position of j in the list P
        for a in range(0,len(self.L)):
            if j in self.L[a]:
                aa=list(self.L[a])
                aa.sort()
                return a,aa.index(j)

        #coloring such that adjacent MS has very different color
        #https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
    def ColorPosition(self,j):#Given a list of coloring P and a vertex j return the color associate to the position of j in the list P
        color_list=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5',
                    '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                    '#cab2d6','#6a3d9a','#ffff99','#b15928']
        jj = list(self.poset.vertices()).index(j)
        return '"' + color_list[jj%20] + '"'

#         #gradient coloring
#         #https://colorbrewer2.org/#type=sequential&scheme=YlGn&n=9
#     def ColorPosition(self,j):#Given a list of coloring P and a vertex j return the color associate to the position of j in the list P
#         (a,aa)=self.PositionInList(j)
#         pink=['#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f']
#         blue=['#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b']
#         green=['#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b']
#         purple=['#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d']
#         red=['#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d']
#         yellow=['#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506']
#         color_list=[green,yellow,blue,purple,pink,red]
#         return '"' + color_list[a%7][aa%8] + '"'


    def __dir__(self):
        return list(self.__dict__.keys()) + dir(self._a)
    def __getattr__(self, attr):
        return getattr(self._a,attr)

    def __init__(self, a, poset, I = set()):
        self._a = a
        self.poset = poset
        self.L = self.poset.ListOfNivelSets()
        self.I=I
        # Compute preimage
        self.preimage_ = {}
        for v in a.complex():
            val = a.value(v)
            if val not in self.preimage_:
                self.preimage_[val] = set()
            self.preimage_[val].add(v)

    def preimage(self, val):
        if val in self.preimage_:
            return self.preimage_[val]
        else:
            return set()

    def graphviz (self):
        """ Return a graphviz string describing the graph and its labels """
        gv = 'digraph {\n'
        indices = { v : str(k) for k,v in enumerate(self.poset.vertices())}
        counts = self._a.count()
        #print(counts)
         ########## minhas mudanças
        n = len(counts[next(iter(counts))])
        def vertex_label(v):
            if v in counts:
                return str(v) + " : " + str(tuple(counts[v]))
            else:
                return str(v) + " : " + str(tuple([0] * n))# return str(v) + " : (0,0,0)"

        for v in self.poset.vertices():

            # self.preimage(v) changed to v in self.poset.vertices()

            if v in self.I:
                gv += indices[v] + '[label="' + vertex_label(v) + ('", shape = septagon, color = red, style=filled, fillcolor=' + self.ColorPosition(v) + '];\n' if v in self.poset.vertices() else '"];\n')



            else:
                gv += indices[v] + '[label="' + vertex_label(v) + ('", style=filled, fillcolor=' + self.ColorPosition(v) + '];\n' if v in self.poset.vertices() else '"];\n')

        if self.L != []: #begin of ranking the nodes that are attractor
            gv += '{rank=same; '

            for v in self.L[0]:
                gv += indices[v] + ' '

            gv += ' ;}; \n' #ending of ranking the nodes that are attractor
            ############

        for v in self.poset.vertices():
            for u in self.poset.children(v):
                gv += indices[v] + ' -> ' + indices[u] + ';\n'

        return gv + '}\n'

    def _repr_svg_(self):
        return graphviz.Source(self.graphviz())._repr_svg_()

