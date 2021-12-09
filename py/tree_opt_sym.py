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
"""
Parameters: 
    G - edge list of original graph; dtype 
    r - root node 
    size - number of tree strategies to generate 
    length - upper bound on tree length

Variables: 
    X - dtype: tupledict, vtype: GRB.BINARY 
        binary variable indicating existance of edge in tree strategy 
    Y - dtype: tupledict, vtype: GRB.BINARY
        binary variable indicating existance of node in tree strategy 
    Z - dtype: tupledict, vtype: GRB.INTEGER, lb: 0, ub = 2^(length-1)
        weight of node i in tree strategy t

Constraints: 
    flow-constr-t-v - num: |T|*|V| 
        sum X[t][(i,j)] = y[t][i] over j in V : (j,i) in A, for all t in T and i in V
        flow constraint to ensure each vertex as in-degree 1
    root-constr-t - num: |T|
        sum Y[t][i] >= 1 over i in V : (r,i) in A, for all t in T
        root constraint to ensure tree starts at root r
    weight-constr-t - num: |V|-1
        sum Z[t][i] >= size over t in T, for all i in V-{r}
        weight constraint to ensure all vertices are included in at leas one tree strategy 
    strat-constr-t-a - |T|*|A|
        Z[t][i] - 2Z[t][j] + 2^length*(1-X[t][(i,j)] >= 0
        weight constraint that guarentees basic strategy, i.e. w(v) = 2w(v+)
    
Objective: 
    Minimize sum of Z[t][i] for all t in T and i in V
"""

# pair tree strategies together in variable definition 
# keeps the same number of constraints (more less?) but does less work since constraints are tighter
def tree_optimization(G, r, size, length):
    # delete spaces in nodes for cartesian product
    # dictionary?
    model = gp.Model('tree-optimizer-symmetry')
    model.setParam(GRB.Param.LogFile, 'lemke_sym-s' + str(size) + '-l' + str(length)+ '.log')
    model.Params.Threads = 1
    print(model.Params.Threads)
    #model = gp.Model('tree-optimizer')

    model._obj = None
    model._bd = None
    model._data = []
    model._start = time.time()

    # model.setParam(GRB.Param.TimeLimit, 4000.0)
    if size%2 != 0: 
        print("Size not divisible by 2. Try again.")
        exit(1)
    size = size//2
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
        model.addConstr(expr - 1 >= 0, name='root-constr-t'+ str(t))

    # Weight Constraint 1
    for i in nodes: # for all i in V\{r}
        if i != r:
            expr = gp.LinExpr()
            for t in range(size): # sum over all t in T
                i_prime = "("+ i[4]+ ", " + i[1] +")"
                expr.addTerms(1.0, Z[t][i])
                expr.addTerms(1.0, Z[t][i_prime])
            model.addConstr(expr - size*2 >= 0, name='weight-constr-t'+ str(t))

    # Strategy Constraint
    for t in range(size): # for all t in T
        for i, j in arcs: # for all (i,j) in A
            if i != r and j != r: # where i,j != r
                #model.addConstr(X[t][(i, j)] * (Z[t][i] - 2*Z[t][j]) == 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))
                model.addConstr(Z[t][i] - 2*Z[t][j] + (2**length)*(1-X[t][(i, j)]) >= 0, name='strat-constr-t'+ str(t)+ '-a'+ str(i)+str(j))

    # Weight Constraint 2 & 3
    for t in range(size): # for all t in T
        for i in nodes: # for all i in V
            model.addConstr(Z[t][i] - (2**(length-1))*Y[t][i] <= 0, name='weight-constr2-t'+ str(t) + '-v' + str(i))
            model.addConstr(Y[t][i] - Z[t][i] <= 0,  name='weight-constr3-t'+ str(t) + '-v' + str(i))

    # # Symmetry Constraints
    # for t in range(size):
    #     for i in nodes: 
    #         #i_prime = (str(i[1]), str(i[0]))
    #         #model.addConstr(X[t][i] - X[t][i_rev] = 0, name='sym-constrX-'+ str(t) + '-v' + str(i))
    #         model.addConstr(X )

    # Set objective
    obj = gp.LinExpr()
    for t in range(size):
        for i in nodes:
            i_prime = "("+ i[4]+ ", " + i[1] +")"
            obj.addTerms(1.0/(size*2), Z[t][i]) 
            obj.addTerms(1.0/(size*2), Z[t][i_prime]) 
    model.setObjective(obj, GRB.MINIMIZE)

    model.optimize(callback=data_cb)

    #model.write('tree_optim' + "_size" + str(size) + "_len" + str(length) + '.lp')
    with open('lemke_sq_sym_log-s60.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(model._data)
    # store results
    trees = {}
    weights = {}
    for t in range(size):
        trees[t] = []
        trees[t+size] = []
        weights[t] = {}
        weights[t+size] = {}
        for i in nodes:
            i_prime = "("+ i[4]+ ", " + i[1] +")"
            weights[t][i] = Z[t][i].x
            weights[t+size][i_prime] = Z[t][i].x
        for i, j in arcs:
            i_prime = "("+ i[4]+ ", " + i[1] +")"
            j_prime = "("+ j[4]+ ", " + j[1] +")"
            if X[t][(i, j)].x:
                trees[t].append((i, j))
                trees[t+size].append((i_prime, j_prime))
        

    save_certificate(weights, trees, n, size, r)
    # Compute pebbling bound
    bound = model.objVal
    bound = math.floor(bound) + 1

    rt = model.Runtime
    print("----------------------------")
    print("|V| = " + str(n))
    print("|E| = " + str(e))
    print("|T| = " + str(size*2))
    print("Max. Length = " + str(length))
    print("Runtime: " + str(rt))
    print("root: " + str(r))
    print("Bound =  " + str(bound))
    print("----------------------------")

    return trees, weights, bound, rt

def cartesian_product(edges1, edges2):
    G = nx.Graph()
    H = nx.Graph()
    G.add_edges_from(edges1)
    H.add_edges_from(edges2)

    prod = nx.cartesian_product(G, H)
    prod_edges = list(prod.edges())

    return prod_edges

def save_certificate(weight_func, paths, n, size, r):


    weights = {}

    for t in paths:
        weights[t] = np.zeros((int(math.sqrt(n)), int(math.sqrt(n))))
        for v in weight_func[t].keys():
            i = int(v[1])
            j = int(v[4])
            weights[t][i][j] = weight_func[t][v]
    for t in weights:
        pd.DataFrame(weights[t]).to_csv("lemke_sq_sym_certificate-t" + str(t) + "-s"+ str(size) + "-v"+ str(r) + ".csv")

    for t in paths:

        pd.DataFrame(paths[t]).to_csv("lemke_sq_sym_edges_tree-t"+  str(t) + "-s"+ str(size) + "-v"+ str(r) + ".csv")

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
    pos = nx.spring_layout(template)

    # plt.subplot(dim_fig, dim_fig, 1)
    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            with_labels=True,
            node_color=color)
    plt.title('Template')
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
                width=list(widths),
                node_size=node_size,
                font_size=font_size,
                labels=weights,
                with_labels=True,
                node_color=color)
        plt.title('Tree Strategy ' + str(t))
        plt.savefig('t' + str(t) + '_size' + str(size) + '_len' + str(length))

    nx.draw(template, pos,
            node_size=node_size,
            font_size=font_size,
            labels=total_weights,
            with_labels=True,
            node_color=color)
    plt.title('Total Weights ')
    plt.savefig('total_size' + str(size) + "_len" + str(length))

def data_cb(model, where):
    if where == gp.GRB.Callback.MIP:
        cur_obj = model.cbGet(gp.GRB.Callback.MIP_OBJBST)
        cur_bd = model.cbGet(gp.GRB.Callback.MIP_OBJBND)

        # Did objective value or best bound change?
        if model._obj != cur_obj or model._bd != cur_bd:
            model._obj = cur_obj
            model._bd = cur_bd
            model._data.append([time.time() - model._start, cur_obj, cur_bd])

def main():
    # Lemke
    # lemke = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
    # r = 6

    # Lemke Square
    lemke = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
    
    lemke_square = cartesian_product(lemke, lemke)
    # visualize_lemke(lemke_square)

    # Bruhat
    # bruhat = [(0,1), (0,2), (0,4), (1,2), (1,7), (2,3), (2,20), (3,23), (4, 5), (4, 8), (5,6), (5,9), (6,7), (6,10), (7,11),
    #          (8,9), (8,16), (9,12), (10,11), (10,13), (11, 19), (12, 13), (12,14), (13,15), (14, 15), (14,17), (15, 18),
    #          (16, 17), (16,20), (17, 21), (18,19), (18,22), (19,23), (20, 21), (21,22), (22,23)]
    # r = 0
    # nontree(bruhat)
    # path = [(0,1), (1,2)]
    # prod_path = cartesian_product(path, path)
    # print(prod_path)
    # r = (0,0)

    # Petersen
    # petersen = [(0,1), (1,2), (2,3), (3,4), (4,0), (0,6), (1,7), (2,8), (3,9), (4,5), (5,7), (5,8), (6,9), (6,8), (7,9)]
    # r = 0

    # Lolipop
    lolipop = [(0,1), (1,2), (2,3), (3,4), (4,1)]
    # r = 0

    test = [((0,0), (1,1)), ((1,1), (2,3)), ((1,1), (3,2)), ((2,3), (4,4)), ((3,2), (4,4))]
    # size = 4
    # length = 6
    # res = tree_optimization (G, r, size, length)
    # trees = res[0]
    # weight_funcs = res[1]
    # G = bruhat
    # sizes = range(6, 10)
    # lengths = range(7, 9)
    # results = {}
    # for s in sizes:
    #     for l in lengths:
    #         res = tree_optimization(G, r, s, l)
    #         trees = res[0]
    #         weight_funcs = res[1]
    #         bound = res[2]
    #         rt = res[3]
    #         results[(s, l)] = [bound, rt]
    #
    #         #visualize_cartesian(G, r, weight_funcs, trees, s, l)
    #
    # save_results(results, sizes, lengths)
    
    G = lemke_square
    # G = [(0,1), (1,2)]
    # G_prod = cartesian_product(G, G)
    s = 20
    l = 10
    r = (0, 0)
    res = tree_optimization(G, r, s, l)
    #trees = res[0]
    #weight_funcs = res[1]
    
    #visualize_cartesian(G_prod, r, weight_funcs, trees, s, l)

    

    # save_certificate(weight_funcs, trees, n)
    # paths_test = [(1,2,1), (1,2)]
    # weight_funcs = {}
    # visualize_cartesian(G, r, weight_funcs, trees, size, length)

if __name__ == "__main__":
    main()