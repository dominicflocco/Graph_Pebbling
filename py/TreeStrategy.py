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

# from PebblingGraph import PebblingGraph 
# from Optimizer import Optimizer 


class TreeStrategy: 

    def __init__(self, graph, root, length): 
        """
        Initializes TreeStrategy Object. 
        
        Attributes: 
            graph - PebblingGraph object 
            root - root of TreeStrategy 
            length - maximum length of TreeStrategy 
            weights - associated weight function that maps vertices to 
                the nonnegative integers 
            nodes - list of nodes in TreeStrategy  
            edges - list of edges in TreeStrategy
        """
        self.root = root 
        self.edges = []
        self.weights = {}
        self.graph = graph 
        self.maxLen = length
        self.nodes = set()
        self.nodes.add(root)
        self.tree = nx.Graph()
        self.diTree = nx.DiGraph()

    def addWeight(self, vertex, weight):
        """
        Adds a weight to the associated weight function of the TreeStrategy 
        
        Parameters: 
            vertex - node in strategy 
            weight - nonnegative integer to be assigned to vertex 
        """
        self.weights[vertex] = weight

    def addEdge(self, src, dst): 
        """ 
        Helper function that adds edges to tree strategy. Used in outputing optimization
        results in Optimization.py.

        Parameters: 
            src - source node of edge 
            dst - destinaation node of edge
        """
        # self.edges.append((str(src), str(dst)))
        # self.nodes.add(str(src))
        # self.nodes.add(str(dst))
         
        self.edges.append((src, dst))
        self.nodes.add(src)
        self.nodes.add(dst)
        # print(src)
    def getWeight(self, vertex):
        """
        Gets the weight of a specified vetex in TreeStrategy 
        
        Parameters: 
            vertex - vertex in TreeStrategy 
        Returns: 
            weight - weight of vertex in TreeStrategy 
        """
        return self.weights[vertex]
    
    def size(self):
        return len(self.edges)

    def saveCertificate(self, filename): 
        """
        For graph products only. Saves certificate ot tree strategy to .csv file 
        in matrix form. Where entry (i, j) represents the weight on vertex (i, j)
        in the tree strategy 

        Parameters: 
            filename - name of .csv file to be written 
        """

        n = int(math.sqrt(len(self.graph.nodes)))
        cert = np.zeros((n,n))
        for v in list(self.nodes):
            i = int(v[2])
            j = int(v[7])
            cert[i][j] = round(self.getWeight(v), 4)   
        pd.DataFrame(cert).to_csv(filename)
    
    def saveEdges(self, filename):
        """
        Saves edges in tree strategy in .csv format for reproducibility. 

        Patameters: 
            filename - name of .csv file to be written 
        """

        pd.DataFrame(self.edges).to_csv(filename)
    
    def visualizeStrategy(self, filename=False):
        """
        Transforms tree stratgy into networkx object and saves a visualization 
        of the strategy to a .png file. 
        
        Parameters: 
            filename - name of .png file to be written. If not specified, 
                        visualization is presented to console but not saved.

        """
        
        self.tree.add_edges_from(self.edges)
        labels = {}
        #print(self.tree.edges)
        # print(self.tree.nodes)
        
        for v in list(self.tree.nodes):
            w = self.getWeight(v)
            label = v + "\n" + str(round(w, 2))
            labels[v] = label
            print(v)
        labels[self.root] = 'r'
        # print(list(self.nodes))
        # Networkx parameters
        node_size = 600
        font_size = 10
        color = 'whitesmoke'
        pos = nx.kamada_kawai_layout(self.tree)

        nx.draw(self.tree, 
                pos=pos, 
                node_size=node_size, 
                font_size=font_size, 
                labels=labels, 
                with_labels=True, 
                node_color=color)
            
        # plt.show()
        if filename: 
            plt.savefig(filename)
        plt.close()

    def verify(self): 
        """ 
        Verifies that the tree strategy outputting by the optimization solver 
        is a valid tree strategy with a valid associated weight function. Does so 
        by checking if the tree is connected and acyclic, and by verifying the weight
        funcion properties are met. 
        """
        valid = True
        self.tree.add_edges_from(self.edges)
        cycles = list(nx.simple_cycles(self.tree))
        if len(cycles) != 0:
            print("Cycle exists in tree. Strategy is invalid.")
            print(self.edges)
            valid = False
        # directedEdges = []
        # visited = set()
        # v = self.root
        # while len(visited) != len(self.nodes):
        #     neighbors = set(self.tree[v])
        #     notvisited = list(neighbors.symmetric_difference(visited))

        #     for u in notvisited:
        #         directedEdges.append((v, u))
        #         visited.add(u)
        return valid
                
    
class NonTreeStrategy: 
    
    def __init__(self, graph, edges):
        """
        Initializes NonTreeStrategy object. 
        
        Attributes: 
            graph - PebblingGraph object 
            root - root of strategy
            weights - associated weight function that maps vertices to 
                the nonnegative integers 
            nodes - list of nodes in TreeStrategy  
            edges - list of edges in TreeStrategy
            lollipop - nx.DiGraph object that stores structure of NonTreeStrategy. 
                Edges are directed towards x0
        
        """
        self.graph = graph 
        self.root = self.graph.root 
        self.edges = list(edges) 
        
        outNodes = set([str(i) for i, j in self.edges])
        inNodes = set([str(j) for i, j in self.edges])
        self.nodes = list(outNodes.union(inNodes))
       
        
        self.lollipop = nx.DiGraph()
  
        self.lollipop.add_edges_from(self.edges)
        
        self.weights = self.generateWeights()
        

    def getWeight(self, vertex):
        """
        Gets the weight of a specified vetex in NonTreeStrategy 
        
        Parameters: 
            vertex - vertex in NonTreeStrategy 
        Returns: 
            weight - weight of vertex in NonTreeStrategy 
        """
        return self.weights[vertex]

    def size(self):
        return len(self.edges)

    def generateWeights(self):
        """
        Generates associated weight function for an even lollipop strategy 
        using Lemma 3.4.1 (Cranston et. al). The function leverages the directed
        edges in even lollipop strategies, which point edges towards x0, to generate
        valid vertex weights. 
        
        """
        # try to make weights integers 

        tail = []

        v = self.graph.root
        while self.lollipop.out_degree(v) == 1:
            u = list(self.lollipop.neighbors(v))[0]
            tail.append((v, u))
            v = u 
        
        xt = v

        
        for v in self.graph.nodes:
            if self.lollipop.in_degree(v) == 2:
                x0 = v 
                break 
      
        v = list(self.lollipop.neighbors(xt))[0]
        u = list(self.lollipop.neighbors(xt))[1]

        path1 = [(xt, v)]
        path2 = [(xt, u)]

        while v != x0: 
            next1 = list(self.lollipop.neighbors(v))[0]
            next2 = list(self.lollipop.neighbors(u))[0]
           
            path1.append((v, next1))
            path2.append((u, next2))
            v = next1
            u = next2
        
        
        t = len(path1)
        s = t + len(tail)
        print(t)
        print(s)
        
        alphaNum = float((2**s + 2**(t-1) - 2))
        alphaDenom = float((2**s - 1))
        alpha = alphaNum/alphaDenom
        weights = {}
        weights[x0] = alpha

        level = 1
        for i in range(1, t+1): 
            v = path1[-i][0]
            u = path2[-i][0]
            weights[v] = 2**level 
            weights[u] = 2**level
            level+=1
        
        tail.reverse()
        for i in range(len(tail)-1): 
            v = tail[i][0]
            weights[v] = 2**level 
            level += 1

        # ensures weights are integers  
        
        if alphaNum%alphaDenom != 0: 
            for key in weights: 
                weights[key] = int(weights[key] * 3)
        
        for v in self.graph.nodes: 
            if v not in weights: 
                weights[v] = 0
    
        return weights

    def trimWeights(self): 
        # use bilevel programing to trim weights
        return None
