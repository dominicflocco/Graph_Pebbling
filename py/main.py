import networkx as nx
import numpy as np
import math
import pandas as pd
import scipy.special
import os 

from Optimizer import Optimizer 
from TreeStrategy import TreeStrategy 
from TreeStrategy import NonTreeStrategy
from PebblingGraph import PebblingGraph 

def countElementaryCirtuits(n):
   
    count = 0 
    for i in range(1, n):
        coef = scipy.special.binom(n, n-i+1)
        count += coef*math.factorial(n-i)
    return count

def treeStrategyOptCrank(graph, size, length, numNodes): 

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "bruhat_crank-cont-v0" + "-s="+ str(size),
                                lpFile = "bruhat_crank-cont-v0" + "-s="+ str(size),
                                paramFile = "bruhat_crank-cont-params" + "-s="+ str(size),
                                objGap=False,
                                timeLimit=False, 
                                cert=False)

    res = optimizer.generateTreeStrategies(graph, size ,length, numNodes)
    
    # save certificates & strategy structure
    strategies = res[0]
    for t in strategies.keys():
        # strategies[t].saveCertificate("bruhat_crank-cont_certv0"+ str(t)+ "-s" + str(size) +".csv")
        strategies[t].saveEdges("bruhat_crank-cont_edgesv0" + str(t)+ "-s" + str(size) + ".csv")

    # visualize strategy 
    for t in strategies.keys():
        strategies[t].visualizeStrategy("bruhat_crank-cont-v0_strategy" + str(t) + "-s" + str(size) + ".png")

def treeStrategyOptSymCrank(graph, size, length, numNodes):

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = None,
                                logFile = "ls_sym2-v" + str(graph.root),
                                lpFile = "ls_sym2-v" + str(graph.root),
                                paramFile = "ls_sym2_params-v" + str(graph.root),
                                objGap=False,
                                timeLimit=False, 
                                cert=True)

    res = optimizer.generateTreeStrategiesSym(graph, size, length, numNodes)
    strategies = res[0]
    for t in strategies.keys():
        strategies[t].saveCertificate("ls_sym2_cert"+ str(t)+ "-v" + str(graph.root) +".csv")
        strategies[t].saveEdges("ls_sym2_edges" + str(t)+ "-v" + str(graph.root) + ".csv")

    # visualize strategy 
    for t in strategies.keys():
        strategies[t].visualizeStrategy("ls_sym2_strategy" + str(t) + "-s" + str(size) + ".png")

def trimStrategies(graph): 

    set_sizes = range(4, 16)
    length = 16

    results = {}
    for size in set_sizes:
        optimizer = Optimizer(solver = "Gurobi", 
                                    threads = 8, 
                                    logFile = "trim_ls-size="+str(size),
                                    lpFile = "trim_ls-size="+str(size),
                                    objGap=False,
                                    timeLimit=43200, 
                                    cert=True)
        
        res = optimizer.generateTreeStrategies(graph, size, length, None)
        strategies = res[0]
        for t in strategies.keys():
            strategies[t].saveCertificate("ls_trim_cert"+ str(t)+ "-s" + str(size) +".csv")
            strategies[t].saveEdges("ls_trim_edges" + str(t)+ "-s" + str(size) + ".csv")
        results[size] = res[1]

    print(results)

def trimStrategiesSym(graph): 

    set_sizes = range(4, 16, 2)
    length = 16

    results = {}
    for size in set_sizes:
        optimizer = Optimizer(solver = "Gurobi", 
                                    threads = 8, 
                                    logFile = "trim_sym_ls-size="+str(size),
                                    lpFile = "trim_sym_ls-size="+str(size),
                                    objGap=False,
                                    timeLimit=43200, 
                                    cert=True)
        
        res = optimizer.generateTreeStrategiesSym(graph, size, length, None)
        strategies = res[0]
        for t in strategies.keys():
            strategies[t].saveCertificate("ls_trim_sym_cert"+ str(t)+ "-s" + str(size) +".csv")
            strategies[t].saveEdges("ls_trim_sym_edges" + str(t)+ "-s" + str(size) + ".csv")
        results[size] = res[1]

    print(results)

def boundLemkeSquare(lemkeSquare): 

    nodes = lemkeSquare.nodes
    bounds = {}
    visited = set()
    # for r in nodes: 
    #     r_rev = "("+ r[4]+ ", " + r[1] +")"
    #     if r_rev not in allNodes and r not in allNodes:
    #         symNodes.append(r)
    #     allNodes.append(r)
    #     allNodes.append(r_rev)
    # rootNodes = symNodes
    # print("Total Number of Nodes: %d" % len(allNodes))
    # print("Number of Symnodes: %d " % len(rootNodes))
    # print(rootNodes)
    params = ["Tree Type", 'Threads', 'Solver', "|V|", '|E|', '|T|', 'Num Nodes', 'Max Len', 'runtime', "root", 'Bound']
    allResults = pd.DataFrame(columns=params)
    #lemkeSquareEdges = list(lemkeSquare.edges())
    size = 12
    length = 16
    for r in nodes:
        if r not in visited: 
            rootdir = "/home/DAVIDSON/doflocco/Graph_Pebbling/results/symmetry/bound-01"
            dir = "/v-" + str(r) 
            dest = rootdir+dir
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            r1, r2 = int(r[1]), int(r[4])
            #lemkeSquare.root = str((r1, r2))
            # print(lemkeSquare.root)
            ls = PebblingGraph(lemkeSquare.edges, (r1, r2))
            optimizer = Optimizer(solver = "Gurobi", 
                                threads = 4, 
                                logFile = "ls_sym_v-" + str(r),
                                lpFile = "ls_sym_v-" + str(r),
                                paramFile = "ls_sym_params" + "-v" + str(r),
                                objGap=False,
                                timeLimit=1800, 
                                cert=True)
            
            res = optimizer.generateTreeStrategiesSym(ls, size, length, None)
            #bounds[r] = res[1]
            strategies = res[0]
            results = res[2] 
            certDest = dest + "/certificates-v" + str(r) 
            os.makedirs(os.path.dirname(certDest), exist_ok=True)
            edgeDest = dest + "/edges-v" + str(r)
            os.makedirs(os.path.dirname(edgeDest), exist_ok=True)

            for t in strategies.keys():
                # strategies[t].saveCertificate(certDest + "/ls_sym_cert"+ str(t)+ "-v" + str(r) +".csv")
                # strategies[t].saveEdges(edgeDest + "/ls_sym_edges" + str(t)+ "-v" + str(r) + ".csv")
                strategies[t].saveCertificate("ls_sym_cert"+ str(t)+ "-v" + str(r) +".csv")
                strategies[t].saveEdges("ls_sym_edges" + str(t)+ "-v" + str(r) + ".csv")
                #strategies[t].visualizeStrategy("ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            # stratDest = dest + "/strategies-v" + str(r)
            # os.makedirs(os.path.dirname(stratDest), exist_ok=True)
            # visualize strategy 
            # for t in strategies.keys():
            #     #strategies[t].visualizeStrategy(stratDest + "/ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            #     strategies[t].visualizeStrategy("ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            allResults = allResults.append(results, ignore_index=True)
            allResults.to_csv(rootdir + "/results.csv")
        visited.add(r)
        r_rev = "("+ r[4]+ ", " + r[1] +")"
        visited.add(r_rev)
    
    # pd.DataFrame.from_dict(bounds, orient="index").to_csv("lemke_square_bounds.csv")


def runNontreeOpt(graph): 

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "lemke_maxUnsolvable.log",
                                lpFile = "lemke_maxUnsolvable.lp",
                                objGap=False,
                                timeLimit = 3600, 
                                cert=True)
    
    lolipops = graph.generateLolipops()  
    strategies = {}
    i = 0
    for pop in lolipops: 
        strategies[i] = NonTreeStrategy(graph, list(pop.edges))
        print(strategies[i].edges)
        print(strategies[i].weights)
        i += 1
    res = optimizer.maxUnsolvable(strategies, graph, graph.root)

    #print(res)

def runHybridOpt(graph):

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "lemke_hybrid.log",
                                lpFile = "lemke_hybrid.lp",
                                objGap=False,
                                timeLimit = False, 
                                cert=True)
    
    lollipops = graph.generateLolipops(cycleLen=4)

    strategies = []
    for pop in lollipops: 
        strategies.append(NonTreeStrategy(graph, list(pop.edges)))
        
    optimizer.hybridOpt(graph, strategies, 5, 5)

def runMaxUnsolvable(graph): 
    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "3cube-maxunsolvable.log",
                                lpFile = "3cube-maxunsolvable.lp",
                                objGap=False,
                                timeLimit = False, 
                                cert=True)
    
    strategies = graph.generateLollipops(cycleLen=4)
    bound = optimizer.maxUnsolvable(strategies, graph)
    print("root = " + str(graph.root) + " -> " + str(bound))

def confirmResults(graph, size, length, numNodes, name): 

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = name + "-s="+ str(size),
                                lpFile = name + "-s="+ str(size),
                                paramFile = name + "-params" + "-s="+ str(size),
                                objGap=False,
                                timeLimit=False, 
                                cert=False)

    res = optimizer.generateTreeStrategies(graph, size, length, numNodes)
    
    # save certificates & strategy structure
    strategies = res[0]
    for t in strategies.keys():
        # strategies[t].saveCertificate("bruhat_crank-cont_certv0"+ str(t)+ "-s" + str(size) +".csv")
        strategies[t].saveEdges(name + "_edges" + str(t)+ "-s" + str(size) + ".csv")

    # # visualize strategy 
    for t in strategies.keys():
        strategies[t].visualizeStrategy(name + "strategy" + str(t) + "-s" + str(size) + ".png")
    #graph.vizualize(strategies, name + "strategy" + str(t) + "-s" + str(size) + ".png", False)

def main(): 
    # Lemke 
    lemke_edges = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
    r = 2
    lemke = PebblingGraph(lemke_edges, r)
    
    #runMaxUnsolvable(lemke)
    lemke_nodes = lemke.nodes
    for r in lemke_nodes:
        lemke = PebblingGraph(lemke_edges, r)
        #runMaxUnsolvable(lemke)
    #strategies = lemke.generateLollipops()
    
    # # Lemke Square 
    lemkeSquare = nx.cartesian_product(lemke.graph, lemke.graph)
    lemkeSquare = list(lemkeSquare.edges())
    r = (0,0)
    #print(countElementaryCirtuits(64))
    lemkeSquare = PebblingGraph(lemkeSquare, r)
    # treeStrategyOpt(lemkeSquare, 6, 15)
    # trimStrategiesSym(lemkeSquare)
    #print(lemkeSquare.nodes)
    #print(lemkeSquare.edges)

    #treeStrategyOptCrank(lemkeSquare, 6, 12, None)
    # boundLemkeSquare(lemkeSquare)
    treeStrategyOptSymCrank(lemkeSquare, 20, 16, None)
    # # Bruhat 
    bruhat = [(0,1), (0,2), (0,4), (1,2), (1,7), (2,3), (2,20), (3,23), (4, 5), (4, 8), (5,6), (5,9), (6,7), (6,10), (7,11),
              (8,9), (8,16), (9,12), (10,11), (10,13), (11, 19), (12, 13), (12,14), (13,15), (14, 15), (14,17), (15, 18),
              (16, 17), (16,20), (17, 21), (18,19), (18,22), (19,23), (20, 21), (21,22), (22,23)]
    r = 1
    # bruhat = PebblingGraph(bruhat, r)
    #treeStrategyOptCrank(bruhat, 8, 8, None)
    # # n-cube
    n = 3
    # hypercube = [(0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)]
    # hypercube = [(0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)]
    # hypercube = [(str(i), str(j)) for i,j in hypercube]
    # r = 1
    # hypercube = PebblingGraph(hypercube, r)

    # strategies = hypercube.generateLollipops(cycleLen=4)
    # for s in strategies:
    #     print(s.weights)
    # runMaxUnsolvable(hypercube)
    # # print(len(hypercube.edges))
    # # print(len(hypercube.nodes))
    # runNontreeOpt(hypercube)

    # Petersen 
    petersen = [(0,1), (1,2), (2,3), (3,4), (4,0), (0,6), (1,7), (2,8), (3,9), (4,5), (5,7), (5,8), (6,9), (6,8), (7,9)]
    r = 5
    # petersen = PebblingGraph(petersen, r)
    # confirmResults(bruhat, 6, 8, None, "bruhat")
    # Lolipop 
    # lolipop = [(0,1), (1,2), (2,3), (3,4), (4,1)]
    # lolipop = PebblingGraph(lolipop, r)

    #lolipops = lemke.generateLolipops(cycleLen=4) 
    # strategies = {}
    # i = 0

    # for pop in lolipops:
    #     #print(pop)
    #     strategies[i] = NonTreeStrategy(lemke, pop.edges)
    #     print(strategies[i].weights)
    #     i+=1
    #boundLemkeSquare(lemkeSquare)

    # optimizer = Optimizer(solver = "Gurobi", 
    #                             threads = 32, 
    #                             logFile = "lemke_square_v-" + str(r),
    #                             lpFile = "lemke_square_v-" + str(r),
    #                             objGap=False,
    #                             timeLimit=False, 
    #                             cert=True)
    
    # res = optimizer.generateTreeStrategiesSym(lemkeSquare, 20, 10)

if __name__ == "__main__":
    main()