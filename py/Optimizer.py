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
import cplex 

from TreeStrategy import TreeStrategy 
from TreeStrategy import NonTreeStrategy
from PebblingGraph import PebblingGraph 


class Optimizer: 

    def __init__(self, solver, logFile, lpFile, threads=False, objGap=False, timeLimit=False, cert=False):
        """
        Initialize the Optimizer class. 

        Parameters: 
            solver - "CPLEX" or "Gurobi" 
            threads - number of threads 
            logFile - string log filename to be saved
            lpFile - string lp filename to be saved
            timeLimit - time limit on optimization solver

        """
        self.solver = solver
        self.threads = threads
        self.timeLimit = timeLimit
        self.lpFile = lpFile
        self.logFile = logFile
        self.objGap = objGap
        self.cert = cert

    def printResults(self, n, e, size, length, rt, r, bound):
        """
        Pretty print results and important information from optimization to console. 

        Parameters: 
            n - number of vertices |V|
            e - number of edges |E|
            size - number of strategies generated
            length - max length of tree strategies generated 
            rt - runtime of solver 
            r - root vertex 
            bound - objective bound 

        """
        
        print("----------------------------")
        print("|V| = " + str(n))
        print("|E| = " + str(e))
        print("|T| = " + str(size*2))
        print("Max. Length = " + str(length))
        print("Runtime: " + str(rt))
        print("root: " + str(r))
        print("Bound =  " + str(bound))
        print("----------------------------")

    def printHybridResults(self, n, e, N, size, length, rt, r, bound):
        """
        Pretty print results and important information from optimization to console. 

        Parameters: 
            n - number of vertices |V|
            e - number of edges |E|
            size - number of strategies generated
            length - max length of tree strategies generated 
            rt - runtime of solver 
            r - root vertex 
            bound - objective bound 

        """
        
        print("----------------------------")
        print("|V| = " + str(n))
        print("|E| = " + str(e))
        print("|T| = " + str(size))
        print("|L| = " + str(N-n))
        print("|S| = " + str(N))
        print("Max. Length = " + str(length))
        print("Runtime: " + str(rt))
        print("root: " + str(r))
        print("Bound =  " + str(bound))
        print("----------------------------")
        
    def hybridOpt(self, graph, lollipops, size, length):

        model = gp.Model('hybrid-optimizer')
        model.setParam(GRB.Param.LogFile, self.logFile + '-s' + str(size) + '-l' + str(length)+ '.log')
        if self.threads:
            model.Params.Threads = self.threads

        if self.timeLimit:
            model.setParam(GRB.Param.TimeLimit, self.timeLimit)

        if self.objGap:
            model.setParam(GRB.Param.MIPGap, self.objGap)

        
        r = str(graph.root)
        n = len(graph.nodes)
        N = n + len(lollipops)
        e = len(graph.edges)
        ubz = 2 ** (length - 1)

        # Add variables to model
        X = {}
        Y = {}
        Z = {}

        arcs = list(graph.arcs)
        nodes = graph.nodes
        # print(arcs)
        # print()
        # print(set(arcs))
        for i in range(size):
            X[i] = model.addVars(arcs, vtype=GRB.BINARY, name='X') 
            Y[i] = model.addVars(nodes, vtype=GRB.BINARY, name='Y')
            Z[i] = model.addVars(nodes, lb=0, ub=ubz, vtype=GRB.INTEGER, name='Z')

        B = model.addVar(range(len(lollipops)), lb=0, vtype=GRB.INTEGER, name="beta")
        # Add Constraints
        # Flow Constraint
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                # sum over all j in V: (j,i) in arcs
                model.addConstr(X[t].sum('*', i) - Y[t][i] == 0, name='flow-constr-t'+str(t) + '-v'+ str(i))

        # Root Constraint
        for t in range(size): # for all t in T
            expr = gp.LinExpr()
            for i, j in graph.arcs:
                if i == r:
                    expr.addTerms(1.0, Y[t][j])
            model.addConstr(expr - 1 >= 0, name='root-constr-t'+ str(t))

        # Weight Constraint 1
        for i in graph.nodes: # for all i in V\{r}
            if i != r:
                expr = gp.LinExpr()
                for t in range(size): # sum over all t in T
                    expr.addTerms(1.0, Z[t][i])
                for l in range(lollipops):
                    expr.addTerms(B[l], l.getWeight(i))
                model.addConstr(expr - N >= 0, name='weight-constr-t'+ str(t) + "-l" + str(l))

        # Strategy Constraint
        for t in range(size): # for all t in T
            for i, j in graph.arcs: # for all (i,j) in A
                if i != r and j != r: # where i,j != r
                    model.addConstr(Z[t][i] - 2*Z[t][j] + (2**length)*(1-X[t][(i, j)]) >= 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))

        # Weight Constraint 2 & 3
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                model.addConstr(Z[t][i] - (2**(length-1))*Y[t][i] <= 0, name='weight-constr2-t'+ str(t) + '-v' + str(i))
                model.addConstr(Y[t][i] - Z[t][i] <= 0,  name='weight-constr3-t'+ str(t) + '-v' + str(i))


        # Set objective
        obj = gp.LinExpr()
        for t in range(size):
            for i in graph.nodes:
                obj.addTerms(1.0/N, Z[t][i]) 
        for l in range(lollipops):
            for i in graph.nodes:
                obj.addTerms(B[i]/N, l.getWeight(i))
        
        model.setObjective(obj, GRB.MINIMIZE)
        lpFilename = self.lpFile + '.lp'
        model.write(lpFilename)

        if self.solver == "CPLEX":
            trees = {} 
            cplexModel = cplex.importModel(lpFilename)
            cplexModel.solve()

            for t in range(size):
                trees[t] = TreeStrategy(graph, graph.root, length)
                trees[t+size] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    trees[t].addWeight(i, cplexModel.GetValue(Z[t][i]))
                    trees[t+size].addWeight(i_prime, cplexModel.GetValue(Z[t][i]))
                for i, j in graph.arcs:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    j_prime = "("+ j[4]+ ", " + j[1] +")"
                    if cplexModel.GetValue(X[t][(i, j)]):
                        trees[t].addEdge(i, j)
                        trees[t+size].addEdge(i_prime, j_prime)
            bound = cplexModel.getObjValue()
            bound = math.floor(bound) + 1

        else: 

            model.optimize()

            # store results
            strategies = {}
            
            for t in range(size):
                strategies[t] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    strategies[t].addWeight(i, Z[t][i].x)
                    
                for i, j in graph.arcs:
                    if X[t][(i, j)].x:
                        strategies[t].addEdge(i, j)

            # save_certificate(weights, trees, n, size, r)
            # Compute pebbling bound
            bound = model.objVal
            bound = math.floor(bound) + 1

            rt = model.Runtime
        
        #self.saveCertificate(strategies)
        self.printHybridResults(n, e, N, size, length, rt, r, bound) 
        
        return strategies, bound, rt
        
    def saveCertificate(self, strategies): 

        weights = {}
        t = 0
        for strat in strategies.keys():
            strategy = strategies[strat]
            weights[t] = np.zeros((int(math.sqrt(len(strategy.nodes))), int(math.sqrt(len(strategy.nodes)))))
            for v in strategy.nodes:
                i = int(v[1])
                j = int(v[4])
                weights[t][i][j] = strategy.weights[v]
            t+= 1

        for t in weights:
            pd.DataFrame(weights[t]).to_csv("lemke_sq_certificate-v" + str(strategies[0].root) + ".csv")
        i = 0

        for t in strategies:
            pd.DataFrame(strategies[t].edges).to_csv("lemke_sq_edges_tree-v"+ str(strategies[0].root) + ".csv")
            i += 1 
    def generateTreeStrategies(self, graph, size, length):
        """
        Generates a set of tree strategies using linear integer programming and a specified
        solver. 

        Paramaters:
            graph - PebblingGraph object 
            size - number of tree strategies to generate
            length - max length of tree 

        Returns: 
            strategies - set of TreeStrategy objectes generated by optimziation 
            bound - pebbling bound found from set of strategies 
            
        """ 
        model = gp.Model('tree-optimizer')

       
        if self.threads: 
            model.Params.Threads = self.threads
        
        if self.timeLimit: 

            model.setParam(GRB.Param.TimeLimit, self.timeLimit)
        
        model.setParam(GRB.Param.LogFile, self.logFile + '-s' + str(size) + '-l' + str(length)+ '.log')
        model.setParam(GRB.Param.MIPGap, self.objGap)
        r = str(graph.root)
        e = len(graph.edges)
        n = len(graph.nodes)

        ubz = 2 ** (length - 1)

        # Add variables to model
        X = {}
        Y = {}
        Z = {}
        
        for i in range(size):
            X[i] = model.addVars(graph.arcs, vtype=GRB.BINARY, name='X') 
            Y[i] = model.addVars(graph.nodes, vtype=GRB.BINARY, name='Y')
            Z[i] = model.addVars(graph.nodes, lb=0, ub=ubz, vtype=GRB.INTEGER, name='Z')

        # Add Constraints
        # Flow Constraint
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                # sum over all j in V: (j,i) in arcs
                model.addConstr(X[t].sum('*', i) - Y[t][i] == 0, name='flow-constr-t'+str(t) + '-v'+ str(i))

        # Root Constraint
        for t in range(size): # for all t in T
            expr = gp.LinExpr()
            for i, j in graph.arcs:
                if i == r:
                    expr.addTerms(1.0, Y[t][j])
            model.addConstr(expr - 1 >= 0, name='root-constr-t'+ str(t))

        # Weight Constraint 1
        for i in graph.nodes: # for all i in V\{r}
            if i != r:
                expr = gp.LinExpr()
                for t in range(size): # sum over all t in T
                    expr.addTerms(1.0, Z[t][i])
                model.addConstr(expr - size >= 0, name='weight-constr-t'+ str(t))

        # Strategy Constraint
        for t in range(size): # for all t in T
            for i, j in graph.arcs: # for all (i,j) in A
                if i != r and j != r: # where i,j != r
                    #model.addConstr(X[t][(i, j)] * (Z[t][i] - 2*Z[t][j]) == 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))
                    model.addConstr(Z[t][i] - 2*Z[t][j] + (2**length)*(1-X[t][(i, j)]) >= 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))

        # Weight Constraint 2 & 3
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                model.addConstr(Z[t][i] - (2**(length-1))*Y[t][i] <= 0, name='weight-constr2-t'+ str(t) + '-v' + str(i))
                model.addConstr(Y[t][i] - Z[t][i] <= 0,  name='weight-constr3-t'+ str(t) + '-v' + str(i))

        # Set objective
        obj = gp.LinExpr()
        for t in range(size):
            for i in graph.nodes:
                obj.addTerms(1.0/size, Z[t][i])

        model.setObjective(obj, GRB.MINIMIZE)

        lpFilename = self.lpFile + "_size" + str(size) + "_len" + str(length) + '.lp'
        model.write(lpFilename)

        if self.solver == "CPLEX": 
            cplexModel = cplex.importModel(lpFilename)
            cplexModel.solve()

            strategies = {}
            for t in range(size):
                trees[t] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    tree[t].addWeight(i, cplexModel.GetValue(Z[t][i]))
                for i, j in graph.arcs:
                    if cplexModel.GetValue(X[t][(i, j)]):
                        trees[t].addEdge(i, j)

            bound = cplexModel.getObjValue()
            bound = math.floor(bound) + 1
        else: 
            model.optimize()
            # store results
            strategies = {}
            for t in range(size):
                trees[t] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    tree[t].addWeight(i, Z[t][i].x)
                for i, j in graph.arcs:
                    if X[t][(i, j)].x:
                        trees[t].addEdge(i, j)

  
            bound = model.objVal
            bound = math.floor(bound) + 1

            rt = model.Runtime
       
        self.printResults(n, e, length, rt, r, bound)

        return strategies, bound
     
    def generateTreeStrategiesSym(self, graph, size, length):
        """
        Generates a set of tree strategies using linear integer programming and a specified
        solver. Leverages symmetry on graph product to reduce constraints. 

        Paramaters:
            graph - PebblingGraph object of cartesian product
            size - number of tree strategies to generate
            length - max length of tree 

        Returns: 
            strategies - set of TreeStrategy objectes generated by optimziation 
            bound - pebbling bound found from set of strategies 
            
        """
        model = gp.Model('tree-optimizer-symmetry')
        model.setParam(GRB.Param.LogFile, self.logFile + '-s' + str(size) + '-l' + str(length)+ '.log')
        if self.threads:
            model.Params.Threads = self.threads

        if self.timeLimit:
            model.setParam(GRB.Param.TimeLimit, self.timeLimit)

        if self.objGap:
            model.setParam(GRB.Param.MIPGap, self.objGap)

        if size%2 != 0: 
            print("Size not divisible by 2. Try again.")
            exit(1)
        size = size//2
        r = str(graph.root)
        n = len(graph.nodes)
        e = len(graph.edges)
        ubz = 2 ** (length - 1)

        # Add variables to model
        X = {}
        Y = {}
        Z = {}

        arcs = list(graph.arcs)
        nodes = graph.nodes
        # print(arcs)
        # print()
        # print(set(arcs))
        for i in range(size):
            X[i] = model.addVars(arcs, vtype=GRB.BINARY, name='X') 
            Y[i] = model.addVars(nodes, vtype=GRB.BINARY, name='Y')
            Z[i] = model.addVars(nodes, lb=0, ub=ubz, vtype=GRB.INTEGER, name='Z')
        
        # Add Constraints
        # Flow Constraint
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                # sum over all j in V: (j,i) in arcs
                model.addConstr(X[t].sum('*', i) - Y[t][i] == 0, name='flow-constr-t'+str(t) + '-v'+ str(i))

        # Root Constraint
        for t in range(size): # for all t in T
            expr = gp.LinExpr()
            for i, j in graph.arcs:
                if i == r:
                    expr.addTerms(1.0, Y[t][j])
            model.addConstr(expr - 1 >= 0, name='root-constr-t'+ str(t))

        # Weight Constraint 1
        for i in graph.nodes: # for all i in V\{r}
            if i != r:
                expr = gp.LinExpr()
                for t in range(size): # sum over all t in T
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    expr.addTerms(1.0, Z[t][i])
                    expr.addTerms(1.0, Z[t][i_prime])
                model.addConstr(expr - size*2 >= 0, name='weight-constr-t'+ str(t))

        # Strategy Constraint
        for t in range(size): # for all t in T
            for i, j in graph.arcs: # for all (i,j) in A
                if i != r and j != r: # where i,j != r
                    model.addConstr(Z[t][i] - 2*Z[t][j] + (2**length)*(1-X[t][(i, j)]) >= 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))

        # Weight Constraint 2 & 3
        for t in range(size): # for all t in T
            for i in graph.nodes: # for all i in V
                model.addConstr(Z[t][i] - (2**(length-1))*Y[t][i] <= 0, name='weight-constr2-t'+ str(t) + '-v' + str(i))
                model.addConstr(Y[t][i] - Z[t][i] <= 0,  name='weight-constr3-t'+ str(t) + '-v' + str(i))


        # Set objective
        obj = gp.LinExpr()
        for t in range(size):
            for i in graph.nodes:
                i_prime = "("+ i[4]+ ", " + i[1] +")"
                obj.addTerms(1.0/(size*2), Z[t][i]) 
                obj.addTerms(1.0/(size*2), Z[t][i_prime]) 
        
        model.setObjective(obj, GRB.MINIMIZE)
        lpFilename = self.lpFile + "_size" + str(size) + "_len" + str(length) + '.lp'
        model.write(lpFilename)

        if self.solver == "CPLEX":
            trees = {} 
            cplexModel = cplex.importModel(lpFilename)
            cplexModel.solve()

            for t in range(size):
                trees[t] = TreeStrategy(graph, graph.root, length)
                trees[t+size] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    trees[t].addWeight(i, cplexModel.GetValue(Z[t][i]))
                    trees[t+size].addWeight(i_prime, cplexModel.GetValue(Z[t][i]))
                for i, j in graph.arcs:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    j_prime = "("+ j[4]+ ", " + j[1] +")"
                    if cplexModel.GetValue(X[t][(i, j)]):
                        trees[t].addEdge(i, j)
                        trees[t+size].addEdge(i_prime, j_prime)
            bound = cplexModel.getObjValue()
            bound = math.floor(bound) + 1

        else: 

            model.optimize()

            # store results
            strategies = {}
            
            for t in range(size):
                strategies[t] = TreeStrategy(graph, graph.root, length)
                strategies[t+size] = TreeStrategy(graph, graph.root, length)
                for i in graph.nodes:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    strategies[t].addWeight(i, Z[t][i].x)
                    strategies[t+size].addWeight(i_prime, Z[t][i].x)
                for i, j in graph.arcs:
                    i_prime = "("+ i[4]+ ", " + i[1] +")"
                    j_prime = "("+ j[4]+ ", " + j[1] +")"
                    if X[t][(i, j)].x:
                        strategies[t].addEdge(i, j)
                        strategies[t+size].addEdge(i_prime, j_prime)

            # save_certificate(weights, trees, n, size, r)
            # Compute pebbling bound
            bound = model.objVal
            bound = math.floor(bound) + 1

            rt = model.Runtime
        
        self.saveCertificate(strategies)
        self.printResults(n, e, size, length, rt, r, bound) 
        
        return strategies, bound, rt
        
    def saveCertificate(self, strategies): 

        weights = {}
        t = 0
        for strat in strategies.keys():
            strategy = strategies[strat]
            weights[t] = np.zeros((int(math.sqrt(len(strategy.nodes))), int(math.sqrt(len(strategy.nodes)))))
            for v in strategy.nodes:
                i = int(v[1])
                j = int(v[4])
                weights[t][i][j] = strategy.weights[v]
            t+= 1

        for t in weights:
            pd.DataFrame(weights[t]).to_csv("lemke_sq_certificate-v" + str(strategies[0].root) + ".csv")
        i = 0

        for t in strategies:
            pd.DataFrame(strategies[t].edges).to_csv("lemke_sq_edges_tree-v"+ str(strategies[0].root) + ".csv")
            i += 1 

    def maxUnsolvable(self, strategies, graph):
        """
        Calculates pebbling bound by finding size of maximum unsolvable configuration. 

        Paramaters: 
            strategies - set of TreeStrategy objects 
            graph - Pebbling Graph object 
            root - root vertex
        Returns: 
            bound - pebbling bound obtained by solver 

        """
        # chceck all vertices are covered with strategies 
        covered = set()
        for s in strategies: 
            for v in s.nodes:
                covered.add(v)
        if len(covered) != len(graph.nodes):
            print("Strategy set is invalid: all vertices must be covered.")
            exit(1)

        model = gp.Model('weight-function-optimizer')
        model.setParam(GRB.Param.LogFile, self.logFile)
        model.Params.Threads = self.threads

        if self.timeLimit:
            model.setParam(GRB.Param.TimeLimit, self.timeLimit)
        
        if self.objGap:
            model.setParam(GRB.Param.MIPGap, self.objGap)
        
        nodes = graph.nodes
        
        
        V = model.addVars(nodes, vtype=GRB.INTEGER, lb=0, ub=2**len(nodes), name='V') 
        #B = model.addVars(strategies, vtype=GRB.CONTINUOUS, lb0=0, name="beta")

        for strat in strategies:
            strategy = strat
            configExpr = gp.LinExpr()
            weightExpr = 0
            for v in nodes: 
                weight = strategy.getWeight(str(v))
                configExpr.addTerms(weight, V[v])
                weightExpr += weight
            model.addConstr(configExpr - weightExpr <= 0)
        
        obj = gp.LinExpr()
            if v != graph.root:

                obj.addTerms(1.0, V[v])
        model.setObjective(obj, GRB.MAXIMIZE)
        model.write(self.lpFile)
        model.optimize()
        
        objVal = model.objVal 
        bound = objVal + 1

        return bound
       

    def generateNonTreeStrategies(self, graph, root, size, length, solver):
        # generate non-tree strategies with optimization 
        return None 