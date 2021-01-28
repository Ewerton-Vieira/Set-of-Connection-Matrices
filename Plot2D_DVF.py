from pychomp import *
import math
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Arc, RegularPolygon
import numpy as np
from numpy import radians as rad
from colour import Color


# Given a Discrete Vector Field, plot Morse Sets and 2D Discrete Vector Field.
class Plot2D_DVF:

    def p_pst(self, angle, radius, x, y):  # function need to draw hexagon
        angle = math.pi * angle / 180
        return (radius * math.cos(angle) + x, radius * math.sin(angle) + y)

    def v_(self, j):  # return dual cells that correspond to fibration value j
        return [i for i in range(0, self.fibration.complex().size()) if self.fibration.value(i) == j]

    # return a collection of cubical cells that correspond to the vertex j in Dag graph
    def DagVtoCell(self, j):
        return [i for i in self.v_(j) if i in self.fibration.complex().topstar(i)]

    def star_dim1(self, k):   # for a given cell return all 1 dimension cells in the star set that has 2 exact topstars
        return [simp for simp in self.simp_dim_1 if len(list(set(self.cc.topstar(k)) & set(self.cc.topstar(simp)))) == 2]

    # draw basic top cell blowup (position, face color, thikness of line, boundary color)
    def TopCellBlowUp(self, ii, jj, color="none", l=1, color2="black", h=' '):
        s_v = 0
        verts = [self.p_pst(i, s_v, ii + k, jj + l) for i, k, l in [(45, 0, 0),
                                                                    (315, 0, 1), (225, 1, 1), (135, 1, 0), (45, 0, 0)]]  # vertices
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        path = Path(verts, codes)
        patch = patches.PathPatch(
            path, facecolor=color, edgecolor=color2, lw=l, hatch=h)
        self.ax.add_patch(patch)


# functions needed for coloring MS

    # Given a list of coloring P and a vertex j return the position of j in the list P

    def PositionInList(self, j):
        for a in range(0, len(self.L)):
            if j in self.L[a]:
                aa = list(self.L[a])
                aa.sort()
                return a, aa.index(j)

        # https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
        # coloring such that adjacent MS has very different color
    # Given a list of coloring P and a vertex j return the color associate to the position of j in the list P
    def ColorPosition(self, j):
        color_list = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
                      '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
                      '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
        jj = list(self.Poset.vertices()).index(j)
        return color_list[jj % 20]

#         #https://colorbrewer2.org/#type=sequential&scheme=YlGn&n=9
#         #gradient coloring
#     def ColorPosition(self,j):#Given a list of coloring P and a vertex j return the color associate to the position of j in the list P
#         (a,aa)=self.PositionInList(j)
#         pink=['#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f']
#         blue=['#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b']
#         green=['#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b']
#         purple=['#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d']
#         red=['#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d']
#         yellow=['#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506']
#         color_list=[green,yellow,blue,purple,pink,red]
#         return color_list[a%7][aa%8]

# end
# #############

    # draw circle arrow
    def drawCirc(self, radius, centX, centY, angle_, theta2_, color_='black'):
        # ========Line
        arc = Arc([centX, centY], radius, radius, angle=angle_,
                  theta1=0, theta2=theta2_, capstyle='round', linestyle='-', lw=2, color=color_)
        self.ax.add_patch(arc)

        # ========Create the arrow head
        # Do trig to determine end position
        endX = centX + (radius / 2) * np.cos(rad(theta2_ + angle_))
        endY = centY + (radius / 2) * np.sin(rad(theta2_ + angle_))

        self.ax.add_patch(  # Create triangle as arrow head
            RegularPolygon(
                (endX, endY),            # (x,y)
                3,                       # number of vertices
                radius / 3,                # radius
                rad(angle_ + theta2_),     # orientation
                color=color_
            )
        )
        # self.ax.set_xlim([centX-radius,centY+radius]) and self.ax.set_ylim([centY-radius,centY+radius])
        # Make sure you keep the axes scaled or else arrow will distort
    #####

    # Color Down sets MS
    def DownCells(self, s):  # return all cells that belong to the downset of the vertex s
        downC = set()
        for a in self.dag.descendants(s):
            downC = downC.union(set(self.DagVtoCell(a)))
        return downC

    def ColorDownMS(self, s):  # color the down sets of given vertex that represent a Morse Set
        corr = "Crimson"  # color
        h = '/'  # hash
        for k in self.DownCells(s):
            [x, y] = self.cc.coordinates(k)
            if self.cc.cell_dim(k) == 0:
                self.PointBlowUp(x, y, "none", 1, corr, h)
            elif self.cc.cell_dim(k) == 2:
                self.TopCellBlowUp(x, y, "none", 1, corr, h)
            elif self.cc.cell_dim(k) == 1 and self.cc.size() / 2 - 1 < k < self.cc.size() - self.cc.size(2):
                self.EdgeVertBlowUp(x, y, "none", 1, corr, h)
            else:
                self.EdgeHorizBlowUp(x, y, "none", 1, corr, h)

    # color interval

    # Given a poset self.Poset induced from self.dag and an convex set I in self.Poset return the correspoding convex set _I_ in self.dag
    def convex_set_in_dag(self, I):
        if not self.Poset.isConvexSet(I):
            return print("I is not a convex set")
        _I_ = set()
        for a in I:
            for b in I:
                if a < b:
                    _I_.update(self.dag.descendants(b).intersection(
                        self.dag.transpose().descendants(a)))  # downset intersection upset
        return _I_

    def IntervalCells(self, I):  # return all cells that belong to the interval I
        IC = set()
        for a in I:
            IC = IC.union(set(self.DagVtoCell(a)))
        return IC

    def ColorMS_I(self, I):  # color  MS(I) all Morse sets in an convex set I in Poset self.Poset
        corr = "Crimson"  # color
        h = '/'  # hash
        _I_ = self.convex_set_in_dag(I)
        for k in self.IntervalCells(_I_):
            [x, y] = std.cc.coordinates(k)
            if self.cc.cell_dim(k) == 0:
                self.PointBlowUp(x, y, "none", 1, corr, h)
            elif self.cc.cell_dim(k) == 2:
                self.TopCellBlowUp(x, y, "none", 1, corr, h)
            elif self.cc.cell_dim(k) == 1 and self.cc.barycenter(k)[0] - 2 * self.cc.coordinates(k)[0] == 1:
                self.EdgeHorizBlowUp(x, y, "none", 1, corr, h)
            else:
                self.EdgeVertBlowUp(x, y, "none", 1, corr, h)

    ##########

    def plotRookField(self):  # ploting the RookField
        for a in range(self.cc.size() - self.cc.size(std.D), self.cc.size() - 1):
            x_temp, y_temp = self.cc.coordinates(a)
            x_temp2, y_temp2 = self.RookField(a)
            if [x_temp2, y_temp2] == [0, 0]:  # plot self edge
                self.drawCirc(.15, x_temp + .5, y_temp + .5,
                              140, 270, color_='m')
                # ax.plot(x_temp+.5,y_temp+.5,color='m',marker=r'$\circlearrowleft$',ms=40)
            if [x_temp2, y_temp2] != [0, 0]:  # plot the vector field
                self.ax.arrow(x_temp + .5 - 0.2 * x_temp2, y_temp + .5 - 0.2 * y_temp2, 0.1 *
                              x_temp2, 0.1 * y_temp2, head_width=0.1, head_length=0.1, fc='m', ec='m')
        self.ax.set_xlim(0, self.sqrt_temp - 1)
        self.ax.set_ylim(0, self.sqrt_temp - 1)

    def FaceGrid(self):
        # simple coloring without considering Morse sets
        for k in self.simp_dim_2:
            if self.FLAG_FRINGE_DINAMICS:

                [x, y] = self.cc.coordinates(k)
                if x < self.domains[0] and y < self.domains[1]:
                    self.TopCellBlowUp(x, y, 'w')

            else:
                if k in self.top_cc_closure():
                    continue
                [x, y] = self.cc.coordinates(k)
                if x < self.domains[0] and y < self.domains[1]:
                    self.TopCellBlowUp(x, y, 'w')

    def coloringMS(self):
        for j in range(len(self.Poset.vertices())):  # coloring Morse sets
            corr = str(self.ColorPosition(list(self.Poset.vertices())[j]))

            for k in self.DagVtoCell(list(self.Poset.vertices())[j]):
                [x, y] = self.cc.coordinates(k)
                self.TopCellBlowUp(x, y, corr)

                # [x, y] = self.cc.coordinates(k)
                # if self.cc.cell_dim(k) == 0:
                #     self.PointBlowUp(x, y, corr)
                # elif self.cc.cell_dim(k) == 2:
                #     self.TopCellBlowUp(x, y, corr)

                # else:
                #     if self.cc.barycenter(k)[0] - 2 * self.cc.coordinates(k)[0] == 1:
                #         self.EdgeHorizBlowUp(x, y, corr)
                #     else:
                #         self.EdgeVertBlowUp(x, y, corr)

    def plotDoubleArrows(self):  # ploting double arrows
        for (iiii, jjjj) in self.digraph.edges():

            if iiii == jjjj:  # ploting self edge
                x_temp, y_temp = self.cc.coordinates(iiii)
                self.ax.plot(x_temp + .5, y_temp + .5, color='red',
                             marker=r'$\circlearrowleft$', ms=20)

            else:
                if self.FLAG_FRINGE_DINAMICS:
                    if (jjjj, iiii) in self.digraph.edges():
                        x_temp, y_temp = self.cc.barycenter(iiii)
                        x_temp, y_temp = x_temp / 2, y_temp / 2
                        x, y = self.cc.barycenter(jjjj)
                        x, y = x / 2, y / 2
                        dx, dy = (x_temp - x) / 4, (y_temp - y) / 4
                        if all([dx * dx < 1 / 15, dy * dy < 1 / 15]):
                            self.ax.annotate('', xy=(x_temp - dx, y_temp - dy), xytext=(x + dx, y + dy),
                                             arrowprops={'arrowstyle': '<->', 'lw': 5, 'ec': 'gray'}, va='center')
                        else:
                            continue
                else:
                    # dont plow double arrow from fringe cells
                    if iiii in self.top_cc_closure() or jjjj in self.top_cc_closure():
                        continue
                    if (jjjj, iiii) in self.digraph.edges():
                        x_temp, y_temp = self.cc.barycenter(iiii)
                        x_temp, y_temp = x_temp / 2, y_temp / 2
                        x, y = self.cc.barycenter(jjjj)
                        x, y = x / 2, y / 2
                        dx, dy = (x_temp - x) / 4, (y_temp - y) / 4
                        if all([dx * dx < 1 / 15, dy * dy < 1 / 15]):
                            self.ax.annotate('', xy=(x_temp - dx, y_temp - dy), xytext=(x + dx, y + dy),
                                             arrowprops={'arrowstyle': '<->', 'lw': 5, 'ec': 'gray'}, va='center')
                        else:
                            continue

    def plotSingleArrows(self):  # ploting arrows that are not double arrow
        for (iiii, jjjj) in self.digraph.edges():
            if not (jjjj, iiii) in self.digraph.edges():  # not print double edges now
                x_temp, y_temp = self.cc.barycenter(iiii)
                x_temp, y_temp = x_temp / 2, y_temp / 2
                x, y = self.cc.barycenter(jjjj)
                x, y = x / 2, y / 2
                dx, dy = (x - x_temp) / 8, (y - y_temp) / 8

                if all([dx * dx < 1 / 63, dy * dy < 1 / 63]):

                    self.ax.arrow(x_temp + 3 * dx, y_temp + 3 * dy, 2 * dx, 2 * dy, head_width=0.1,
                                  head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                else:
                    if self.FLAG_FRINGE_DINAMICS:
                        self.ax.arrow(x_temp, y_temp, dx / 4, dy / 4, head_width=0.1,
                                      head_length=0.1, fc='g', ec='g', overhang=self.ohang)

                    else:
                        if x_temp == x:
                            self.ax.arrow(x_temp + 3 * dx, y_temp - self.domains[1] + 3 / 8, 2 * dx, 2 / 8, head_width=0.1,
                                          head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                            self.ax.arrow(x_temp + 3 * dx, y_temp - 3 / 8, -2 * dx, -2 / 8, head_width=0.1,
                                          head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                        else:
                            if iiii - jjjj > self.domains[0] * self.domains[1] - 2:
                                self.ax.arrow(x_temp - self.domains[0] + 3 / 8, y_temp - self.domains[1] + 1, 2 / 8, 0, head_width=0.1,
                                              head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                            elif jjjj - iiii > self.domains[0] * self.domains[1] - 2:
                                self.ax.arrow(x_temp - 3/8, y_temp, -2 / 8, 0, head_width=0.1,
                                              head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                            else:
                                self.ax.arrow(x_temp - self.domains[0] + 3 / 8, y_temp + 1, 2 / 8, 0, head_width=0.1,
                                              head_length=0.1, fc='k', ec='k', overhang=self.ohang)

                                self.ax.arrow(x - self.domains[0] + 1 - 3 / 8, y + 1, - 2 / 8, 0, head_width=0.1,
                                              head_length=0.1, fc='k', ec='k', overhang=self.ohang)  # for self.Att

    def __init__(self, std, Poset=None, I=set(), Show_Fringe=False):

        self.FLAG_FRINGE_DINAMICS = Show_Fringe  # True to show the fringe cells

        self.D = std.D  # size
        self.domains = std.cc.boxes()  # a list of numbers of thresholds

        self.cc = std.cc  # Cubical Complex

        self.top_cc_closure = std.top_cc_closure

        self.digraph = std.digraph  # Directed Acyclic Graph()
        self.blowupfringe = std.blowupfringe
        self.blowupinfinity = std.blowupinfinity
        self.Poset = Poset  # poset
        self.I = I  # interval

        (self.dag, self.fibration) = FlowGradedComplex(
            std.complex(), std.diagram())

        self.fig, self.ax = plt.subplots(
            figsize=(3 * self.domains[0], 3 * self.domains[1]))

        self.sqrt_temp = math.sqrt(self.cc.size(0))

        self.ohang = .5  # overhang for the single arrows

        self.simp_dim_1 = [simp for simp in self.cc if self.cc.cell_dim(
            simp) == 1]  # all 1 dim simp

        self.simp_dim_2 = [simp for simp in self.cc if self.cc.cell_dim(
            simp) == 2]  # all 2 dim simp

        self.FaceGrid()

        if Poset:
            self.L = self.Poset.ListOfNivelSets()
            self.coloringMS()

        self.plotDoubleArrows()

        self.plotSingleArrows()

        if self.FLAG_FRINGE_DINAMICS:
            self.ax.set_xlim(0, self.domains[0])
            self.ax.set_ylim(0, self.domains[1])

        else:
            self.ax.set_xlim(-0.5, self.domains[0] - 0.5)
            self.ax.set_ylim(-0.5, self.domains[1] - 0.5)
            plt.axis('off')

        # Morse set of an interval
        if set(I) != set():
            if not self.Poset.isConvexSet(I):
                return print("===== I is not a convex set =====")
            else:
                self.ColorMS_I(self.I)
