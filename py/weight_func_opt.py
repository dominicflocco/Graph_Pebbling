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

def generateNontreeOpt(G, r, size, length): 
    model = gp.Model('nontree-generator')
    model.setParam(GRB.Param.LogFile, 'nontree-generator.log')

    r = str(r)
    e = len(G)

    # create bijection
    arcs = []
    for u, v in G:
        arcs.append((str(v), str(u)))
        arcs.append((str(u), str(v)))
    nodes = list(set([i for i, j in arcs]))
    n = len(nodes)
    ubz = 2 ** (length - 1)

    # Add variables to model
    X = {}
    Y = {}
    Z = {}
    for i in range(size):
        X[i] = model.addVars(arcs, vtype=GRB.BINARY, name='X') 
        Y[i] = model.addVars(nodes, vtype=GRB.BINARY, name='Y')
        Z[i] = model.addVars(nodes, lb=0, ub=ubz, vtype=GRB.INTEGER, name='Z')

    # Add Constraints
    # Flow Constraint
    for t in range(size): # for all t in T
        for i in nodes: # for all i in V
            # sum over all j in V: (j,i) in arcs
            model.addConstr(X[t].sum('*', i) - Y[t][i] == 0, name='flow-constr-t'+str(t) + '-v'+ str(i))

    # Root Constraint
    for t in range(size): # for all t in T
        expr = gp.LinExpr()
        for i, j in arcs:
            if i == r:
                expr.addTerms(1.0, Y[t][j])
        model.addConstr(expr - 1 == 0, name='root-constr-t'+ str(t))

    # Weight Constraint 1
    for i in nodes: # for all i in V\{r}
        if i != r:
            expr = gp.LinExpr()
            for t in range(size): # sum over all t in T
                expr.addTerms(1.0, Z[t][i])
            model.addConstr(expr - size >= 0, name='weight-constr-t'+ str(t))

    # Strategy Constraint
    for t in range(size): # for all t in T
        for i, j in arcs: # for all (i,j) in A
            if i != r and j != r: # where i,j != r
                #model.addConstr(X[t][(i, j)] * (Z[t][i] - 2*Z[t][j]) == 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))
                model.addConstr(Z[t][i] - 2*Z[t][j] + (2**length)*(1-X[t][(i, j)]) >= 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))

    # Quasi-tree Constraint
    for t in range(size): 
        expr1 = gp.LinExpr()
        expr2 = gp.LinExpr()
        for i,j in arcs: 
            if i != r and j != r:
                expr1.addTerms(1.0, X[t][(i,j)])
                expr1.addTerms(1.0, X[t][(j,i)])
                expr2.addTerms(1.0, X[t][(i,j)])
                expr2.addTerms(1.0, Y[t][i])
                expr2.addTerms(-1.0, Y[t][j]) 

        model.addConstr(expr1 - 2*Y[t].sum() + 1 == 0)
        model.addConstr(expr2 - X[t].sum() + 1 == 0)


    # Weight Constraint 2 & 3
    for t in range(size): # for all t in T
        for i in nodes: # for all i in V
            model.addConstr(Z[t][i] - (2**(length-1))*Y[t][i] <= 0, name='weight-constr2-t'+ str(t) + '-v' + str(i))
            model.addConstr(Y[t][i] - Z[t][i] <= 0,  name='weight-constr3-t'+ str(t) + '-v' + str(i))

    # Set objective
    obj = gp.LinExpr()
    for t in range(size):
        for i in nodes:
            obj.addTerms(1.0/size, Z[t][i])

    model.setObjective(obj, GRB.MINIMIZE)

    model.optimize()

    model.write('nontree-generator' + "_size" + str(size) + "_len" + str(length) + '.lp')

    # store results
    trees = {}
    weights = {}
    for t in range(size):
        trees[t] = []
        weights[t] = {}
        for i in nodes:
            weights[t][i] = Z[t][i].x
        for i, j in arcs:
            if X[t][(i, j)].x:
                trees[t].append((i, j))

    # save_certificate(weights, trees, n, size, r)
    # Compute pebbling bound
    bound = model.objVal
    bound = math.floor(bound) + 1

    rt = model.Runtime
    print("----------------------------")
    print("|V| = " + str(n))
    print("|E| = " + str(e))
    print("|T| = " + str(size))
    print("Max. Length = " + str(length))
    print("Runtime: " + str(rt))
    print("root: " + str(r))
    print("Bound =  " + str(bound))
    print("----------------------------")

    return trees, weights, bound, rt

def maxUnsolvableOpt(strategies, weightFunctions, G, r): 
    model = gp.Model('weight-function-optimizer')
    model.setParam(GRB.Param.LogFile, 'weight_optimizer.log')

    nodes = G.nodes

    V = {}
    sum = gp.LinExp()
    
    V = model.addVars(nodes, vtype=GRB.INTEGER, name='V') 
    
    for i in strategies:
        expr = gp.LinExpr()
        weightSum = 0
        for v in nodes: 
            expr.addTerms(weightFunctions[i][v], V[v])
            weightSum += weightFunctions[i][v]
        model.addConstr(expr - weightSum <= 0)

    sum = gp.LinExp()
    for v in nodes: 
        if v != r:
            sum.addTerms(1.0, V[v])
    model.setObjective(sum, GRB.MAXIMIZE)

    model.optimize()
    
    objVal = model.objVal 
    bound = objVal + 1

    return bound

def generateStrategies(G, r): 
    orig = nx.Graph()
    orig.add_edges_from(G)
    
    arcs = [(v,u) for u,v in orig.edges]
    arcs = arcs + list(orig.edges)
    biDirect = nx.DiGraph()
    biDirect.add_edges_from(arcs)

    cycles = nx.simple_cycles(biDirect)
    for c in cycles: 
        if r in c: 
            cycles.remove(c)

    strategies = [shortestPath(G, c, r) for c in cycles]

    weightFunctions = {}
    strategiesDict = {}
    k = 2
    i = 0
    # get length of cycle and divide by 2 
    # add this legnth/2 to the tail length 
    # first weight = k**length 
    # traverse cycle by levels with weight = k**level 
    for strat in strategies: 
        S = nx.DiGraph()
        S.add_edges_from(strat[0])
        lenCycle = strat[1]
        lenTail = strat[2]
        level = lenCycle//2 + lenTail  
        v = r
        weightFunc = {}
        while len(list(S.out_edges(v))) <= 1:
            weightFunc[v] = 2**level
            level =- 1
            v = list(S.out_edges(v))[0]
        splitVertex = v
        left = list(S.out_edges(splitVertex))[0]
        right = list(S.out_edges(splitVertex))[1]
        for i in range(lenCycle//2): 
            weightFunc[left] = 2**level
            weightFunc[right] = 2**level
            left = list(S.out_edges(splitVertex))[0]
            right = list(S.out_edges(splitVertex))[1]
            level =- 1
        weightFunctions[i] = weightFunc
        strategiesDict[i] = strat[0]
        i += 1

    return strategies 
    

def shortestPath(G, cycle, r): 
    
    minLen = float('inf')
    shortestPath = []
    for v in cycle: 
        p = nx.shortest_path(G, v, r)
        if len(p) < minLen:
            shortestPath = p

    tailEdges = [(shortestPath[i], shortestPath[i+1]) for i in range(len(shortestPath) - 1)]
    
    cycleEdges = [(cycle[i], cycle[i+1]) for i in range(len(cycle) -1)]
    
    strategy = tailEdges + cycleEdges
    strategyTup = (strategy, len(cycle), len(tailEdges))
    return strategyTup 
    
# def generateTrees(): 

def visualize(G, r, weight_funcs, paths, size, length):

    # Save copy of template graph
    G = [(str(i), str(j)) for i,j in G]
    r = str(r)
    template = nx.Graph()
    template.add_edges_from(G)
    template = nx.relabel_nodes(template, {r: 'r'})

    # Networkx parameters
    node_size = 100
    font_size = 7
    color = 'whitesmoke'
    dim_fig = math.ceil(math.sqrt(len(paths)))+1
    pos = nx.kamada_kawai_layout(template)

    plt.subplot(dim_fig, dim_fig, 1)
    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            with_labels=True,
            node_color=color)
    plt.title('Template')

    total_weights = {i: 0 for i in template.nodes}
    for t in range(len(paths)):
        orig = nx.Graph()

        orig.add_edges_from(template.edges)

        tree = nx.Graph()
        tree.add_edges_from(paths[t])
        tree = nx.relabel_nodes(tree, {r: 'r'})

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
            elif v not in weight_funcs[t]:
                weight[v] = ""
            else:
                
                weight[v] = int(math.ceil(weight_funcs[t][v]))
                total_weights[v] += int(math.ceil(weight[v]))
                
        nx.set_node_attributes(orig, weight, 'weight')

        weights = nx.get_node_attributes(orig, 'weight')
        
        widths = [orig[u][v]['width'] for u, v in orig.edges()]
        
        plt.subplot(dim_fig, dim_fig, t+2)

        nx.draw(orig, pos,
                width=widths,
                node_size=node_size,
                font_size=font_size,
                labels=weights,
                with_labels=True,
                node_color=color)
        plt.title('Tree Strategy ' + str(t))



    plt.subplot(dim_fig, dim_fig, dim_fig**2)
    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            labels=total_weights,
            with_labels=True,
            node_color=color)
    plt.title('Total Weights')
    plt.savefig('lemke-nontree-generator')

    plt.show()
    plt.close()

def visualize_cartesian(G, r, weight_funcs, paths, size, length):

    # Save copy of template graph
    G = [(str(i), str(j)) for i, j in G]
    r = str(r)
    template = nx.Graph()
    template.add_edges_from(G)
    template = nx.relabel_nodes(template, {r: 'r'})

    # Networkx parameters
    node_size = 200
    font_size = 8
    color = 'whitesmoke'
    # dim_fig = math.ceil(math.sqrt(len(paths)))+1
    pos = nx.kamada_kawai_layout(template)

    # plt.subplot(dim_fig, dim_fig, 1)
    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            with_labels=True,
            node_color=color)
    plt.title('Template')
    plt.show()
    plt.savefig('template.png', format='PNG')

    total_weights = {i: 0 for i in template.nodes}
    for t in range(len(paths)):
        orig = nx.Graph()

        orig.add_edges_from(template.edges)

        tree = nx.Graph()
        tree.add_edges_from(paths[t])
        tree = nx.relabel_nodes(tree, {r: 'r'})

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
            elif v not in weight_funcs[t]:
                weight[v] = ""
            else:
                weight[v] = math.ceil(weight_funcs[t][v])
                total_weights[v] += math.ceil(weight[v])

        nx.set_node_attributes(orig, weight, 'weight')

        weights = nx.get_node_attributes(orig, 'weight')
        
        widths = [orig[u][v]['width'] for u, v in orig.edges()]
        nx.draw(orig, pos,
                width=widths,
                node_size=node_size,
                font_size=font_size,
                labels=weights,
                with_labels=True,
                node_color=color)
        plt.title('Tree Strategy ' + str(t))
        plt.show()
        plt.savefig('t' + str(t) + '_size' + str(size) + '_len' + str(length))

    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            labels=total_weights,
            with_labels=True,
            node_color=color)
    plt.title('Total Weights ')
    plt.show()
    plt.savefig('total_size' + str(size) + "_len" + str(length))
    plt.close()

def main(): 
    # Lemke
    lemke = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
    
    r = 6
    s = 5
    l = 4
    res = generateNontreeOpt(lemke, r, s, l)
    trees = res[0]
    weight_funcs = res[1]
    # print(trees)
    visualize(lemke, r, weight_funcs, trees, s, l)


if __name__ == "__main__":
    main()