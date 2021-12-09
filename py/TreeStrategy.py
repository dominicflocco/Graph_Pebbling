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
        self.root = root 
        self.edges = []
        self.weights = {}
        self.graph = graph 
        self.maxLen = length
        self.nodes = list(set([i for i, j in self.edges]))

    def addWeight(self, vertex, weight): 
        self.weights[vertex] = weight
    
    def addEdge(self, src, dst): 
        self.edges.append((src, dst))

    def getWeight(self, vertex):
        return self.weights[vertex]
    
    def size(self):
        return len(self.edges)

class NonTreeStrategy: 
    
    def __init__(self, graph, edges): 
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
        return self.weights[vertex]

    def size(self):
        return len(self.edges)

    def generateWeights(self):
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
