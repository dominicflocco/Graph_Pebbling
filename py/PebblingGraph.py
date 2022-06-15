import sys
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
import time

# from Optimizer import Optimizer 
from TreeStrategy import TreeStrategy 
from TreeStrategy import NonTreeStrategy

class PebblingGraph:

    def __init__(self, edges, root):
        """
        Initialize PebblingGraph object.

        Parameters:
            edges - list of undirected edges in graph structure
            root - specified root vertex

        Attributes:
            root - root vertex
            graph - networkx undirected graph object
            arcs - bidirected edges
            nodes - list of nodes in graph
            edges - list of undirected edges in graph
            directed - networkx directed graph object with bidirected edges

        """

        self.root = str(root)
        #self.edges = [(str(u), str(v)) for u, v in edges]
        self.edges = list(set(edges))

        self.graph = nx.Graph()
        self.graph.add_edges_from(self.edges)

        self.arcs = []
        for u, v in edges:
            self.arcs.append((str(v), str(u)))
            if u != v:
                self.arcs.append((str(u), str(v)))
        
        self.arcs = list(set(self.arcs))
        
        self.nodes = list(set([str(i) for i, j in self.arcs]))

        self.directed = nx.DiGraph()
        self.directed.add_edges_from(self.arcs)

    def cartesianProduct(self, edges2):
        """
        Generates cartesian product of PebblingGraph object and set of edges using
        built-in networkx cartesian product function.

        Parameters:
            edges2 - list of edges in graph H

        Returns:
            prodEdges - list of edges in resulting cartesian product
        """
        G = self.graph
        H = nx.Graph()
        H.add_edges_from(edges2)

        prod = nx.cartesian_product(G, H)
        prodEdges = list(prod.edges())

        return prodEdges

    def generateAllCycles(self):
        """
        Generates all simple cycles of even length in PebblingGraph object that
        do not contain the root.

        Returns:
            cycles - list of cycles found in PebblingGraph, where each element is
            a list of ordered vertices
        """
        # Add paramater to fix number of cycles generated
        candCycles = list(nx.simple_cycles(self.directed))
        cycles = []
        allCycleEdges = []
        for c in candCycles:
           
            cEdges = [(c[i], c[i+1]) for i in range(len(c) - 1)]
            cEdges.append((c[-1], c[0]))

            if ((self.root not in c) and (len(cEdges)%2 == 0) and 
                (len(cEdges) > 2 ) and (self.isUnique(cEdges, allCycleEdges))):
                cycles.append(c)
                allCycleEdges.append(cEdges)


        return cycles

    def generateCycles(self, length):
        """
        Generates all simple cycles of given length in PebblingGraph objected that do
        not contain the root.

        Returns:
            cycles - list of cycles found in PebblingGraph, where each element is
            a list of ordered vertices
        """
        assert length%2==0, "cycles must have even length"

        candCycles = list(nx.simple_cycles(self.directed))
        cycles = []
        allCycleEdges = []
        for c in candCycles:
           
            cEdges = [(c[i], c[i+1]) for i in range(len(c) - 1)]
            cEdges.append((c[-1], c[0]))
            
            if (self.root not in c) and (len(cEdges) == length) and self.isUnique(cEdges, allCycleEdges):
                
                cycles.append(c)
                allCycleEdges.append(cEdges)

                
        return cycles

    def getTail(self, cycle):
        """
        Generates tail for a given cycle by finding the shortest path from one vertex
        to the root such that the path is disjoint from the cycle.

        Parameters:
            cycle - list of vertices in cycle

        Returns:
            tail - list of edges in tail
        """
        minLen = float('inf')
        shortestPath = []
        
        cycleEdges = [(cycle[i], cycle[i+1]) for i in range(len(cycle) - 1)]
        cycleEdges.append((cycle[-1], cycle[0]))
        #print(cycleEdges)
        for v in cycle:
            
            path = nx.shortest_path(self.graph, v, self.root)
            pathEdges = [(path[i], path[i+1]) for i in range(len(path) - 1)]
            vertices = list(path)
           
            vertices.remove(v)
            
            valid = True
            for e in cycleEdges:
                e_rev = (e[1], e[0])
                if e in pathEdges or e_rev in pathEdges or (e[1] in vertices) or (e[0] in vertices):
                    valid = False
                    break
            
            if len(pathEdges) < minLen and valid:
                shortestPath = pathEdges
                minLen = len(shortestPath)
      
        tailEdges = shortestPath
        directedEdges = []
        for e in tailEdges: 
            directedEdges.append((e[1], e[0])) 
            
        tail = nx.DiGraph() 
        tail.add_edges_from(directedEdges)

        return tail

    def generateLollipops(self, cycleLen=False):
        """
        Constructs lollipop strategies on PebblingGraph object by combining cycles
        and tails.

        Parameters:
            cycleLen - if an integer value is specific, generates all lollipops with given
            cycle length. Otherwise, produces all lollipops.

        Returns:
            lollipops - list of lollipops, where each lollipop is a list of edges
        """
        if not cycleLen:
            cycles = self.generateAllCycles()
        else:
            cycles = self.generateCycles(cycleLen)

        lollipops = []

        for c in cycles:
            tail = self.getTail(c)
            
            for v in tail.nodes: 
                if tail.out_degree(v) == 0: 
                    xt = v
           
            cycle = [(c[i], c[i+1]) for i in range(len(c) - 1)]
            cycle.append((c[-1], c[0]))
           
            cycleEdges = []
            v = xt 
        
            lenCycle = len(cycle)
            for i in range(lenCycle):
                if i == lenCycle/2:
                    v = xt
                for e in cycle: 
                    if e[0] == v: 
                        cycleEdges.append((e[0], e[1]))
                        cycle.remove(e)
                        v = e[1]
                        break
                    elif e[1] == v: 
                        cycleEdges.append((e[1], e[0]))
                        cycle.remove(e)
                        v = e[0] 
                        break


            directedEdges = list(tail.edges) + cycleEdges
            print(list(tail.edges))
            print(cycleEdges)
            strategy = NonTreeStrategy(self, directedEdges)
            lollipops.append(strategy)
      
            
        return lollipops

    def isUnique(self, cycle, allCycles): 
        """
        Determines whethere given cycle is unique from a set of cycles previously
        found in a Pebbling Graph. A cycle is unique from a set of cycles if 
        if has vertex and edge sets that are unique from all previously found cycles. 
        
        Parameters: 
            cycle - object we wish to determine if unique 
            allCycles - set of cycles already found in graph 
        Returns: 
            True if cycle is unique from set all Cycles 
            False otherwise
        """

        uniqueList = []
        arcs = list(set([(j,i) for i,j in cycle]))
        if len(allCycles) == 0: 
            return True
        for cyc in allCycles: 
            unique = False
            candArcs = list(set([(j,i) for i,j in cyc]))
            # if len(arcs) != len(candArcs):
            #     unique = True
            #     uniqueList.append(unique)
            #     continue
            # if len(set(candArcs).intersection(set(arcs))) != len(arcs):
            #     unique = True
            for e in cyc: 
                e_rev = (e[1], e[0])
                if (e not in cycle) and (e_rev not in cycle): 
                    unique = True 
            uniqueList.append(unique)
        
        unique = True
        for u in uniqueList: 
            if not u: 
                unique = False

        return unique 

    def generateSubgraphLollipops(self, strategy):
        """
        Generates lollipop strategies by finding cycles in the subgraph induced by
        a tree strategy.

        Parameters:
            strategy - TreeStrategy object

        Returns:
            lollipops - list of lollipops, where each lollipop is a list of edges
        """

        subgraph = self.directed.subgraph(strategy.nodes).copy()

        cycles = nx.simple_cycles(subgraph)
        for c in cycles:
            cEdges = set([(c[i], c[i+1]) for i in range(len(c) - 1)])
            if (self.root in c) or (len(cEdges)%2 != 0):
                cycles.remove(c)

        lollipops = []
        for c in cycles:
            minLen = float('inf')
            shortestPath = []
            cycleEdges = set([(c[i], c[i+1]) for i in range(len(c) - 1)])

            for v in cycle:
                path = nx.shortest_path(self.graph, v, self.root)
                
                pathEdges = set([(path[i], path[i+1]) for i in range(len(path) - 1)])
                valid = True
                for e in cycleEdges:
                    if (e in pathEdges) or (e[1] in list(path)) or (e[0] in list(path)):
                        valid = False

                if len(pathEdges) < minLen and valid:
                    shortestPath = pathEdges

            tail = shortestPath
            cycle = [(c[i], c[i+1]) for i in range(len(c) - 1)]
            lollipop = tail + cycle
            lollipops.append(lollipop)

        return lollipops

    def vizualize(self, strategies, filename, multiview=True):

        self.graph = nx.relabel_nodes(self.graph, {self.root: 'r'})

        # Networkx parameters
        node_size = 100
        font_size = 7
        color = 'whitesmoke'

        if multiview:
            dim_fig = math.ceil(math.sqrt(len(strategies)))+1
            plt.subplot(dim_fig, dim_fig, 1)

        pos = nx.kamada_kawai_layout(self.graph)

        nx.draw(self.graph, pos,
                node_size=node_size,
                font_size=font_size,
                with_labels=True,
                node_color=color)
        plt.title('Template')

        if not multiview:
            plt.savefig('template.png', format='PNG')

        total_weights = {i: 0 for i in self.graph.nodes}
        for t in range(len(strategies)):
            orig = nx.Graph()

            orig.add_edges_from(self.graph.edges)

            tree = nx.Graph()
            tree.add_edges_from(strategies[t].edges)
            tree = nx.relabel_nodes(tree, {self.root: 'r'})

            width = {}
            for u, v in orig.edges():
                if tree.has_edge(u, v) or tree.has_edge(v,u):
                    width[(u, v)] = 3.0
                    width[(v, u)] = 3.0
                else:
                    width[(u, v)] = 1.0
                    width[(v, u)] = 1.0

            nx.set_edge_attributes(orig, width, 'width')
            weight = {}
            for v in list(orig.nodes):
                if v == 'r':
                    weight['r'] = 'r'
                    total_weights[v] = 'r'
                elif v not in strategies[t].weights:
                    weight[v] = ""
                else:
                    weight[v] = int(math.ceil(strategies[t].getWeight(v)))
                    total_weights[v] += int(math.ceil(weight[v]))

            nx.set_node_attributes(orig, weight, 'weight')

            weights = nx.get_node_attributes(orig, 'weight')

            widths = [orig[u][v]['width'] for u, v in orig.edges()]
            if multiview:
                plt.subplot(dim_fig, dim_fig, t+2)

            nx.draw(orig, pos,
                    width=widths,
                    node_size=node_size,
                    font_size=font_size,
                    labels=weights,
                    with_labels=True,
                    node_color=color)
            plt.title('Tree Strategy ' + str(t))
            if not multiview:
                plt.show()
                plt.savefig('t' + str(t) + '_size' + str(len(strategies)) + '_len' + str(strategies[0].maxLen))
                plt.close()

        #plt.subplot(dim_fig, dim_fig, dim_fig**2)
        nx.draw(self.graph, pos,
                node_size=node_size,
                font_size=font_size,
                labels=total_weights,
                with_labels=True,
                node_color=color)
        plt.title('Total Weights')
        plt.show()
        if multiview:
            plt.savefig(filename)
        else:
            plt.savefig('total_size' + str(len(strategies)) + "_len" + str(strategies[0].maxLen))

        plt.close()

        return None

    def calculateBound(self, strategies):
        """ 
        Calculates the pebbling bound from a given set of strategies manually by 
        finding the average vertex weight across strategies (Covering Lemma). 

        Parameters: 
            strategies - set of TreeStrategy Objects
        """

        sumWeights = 0
        for t in range(len(strategies)):
            strategy = strategies[t]
            for i in self.nodes:
                sumWeights += strategy.getWeight(i)

        avg = sumWeights/len(strategies)
        bound = math.floor(avg) + 1

        return bound
