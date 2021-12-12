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
        self.nodes = list(set([i for i, j in self.edges]))

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
        Adds a directed edge to the TreeStrategy structure. Directed edges 
        point away from root r.
        
        Parameters: 
            src - source vertex of edge to be added
            dst - destination vertex of edge to be added
        """
        self.edges.append((src, dst))

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
