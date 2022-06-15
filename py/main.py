import networkx as nx
import numpy as np
import math
import pandas as pd
import scipy.special
import collections
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
                                logFile = "ls_sym-v" + str(graph.root),
                                lpFile = "ls_sym-v" + str(graph.root),
                                paramFile = "ls_sym_params-v" + str(graph.root),
                                objGap=False,
                                timeLimit=False, 
                                cert=True)

    res = optimizer.generateTreeStrategiesSym(graph, size, length, numNodes)
    strategies = res[0]
    for t in strategies.keys():
        strategies[t].saveCertificate("ls_sym_cert"+ str(t)+ "-v" + str(graph.root) +".csv")
        strategies[t].saveEdges("ls_sym_edges" + str(t)+ "-v" + str(graph.root) + ".csv")

    # visualize strategy 
    # for t in strategies.keys():
    #     strategies[t].visualizeStrategy("ls_sym_strategy" + str(t) + "-s" + str(size) + ".png")

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
    visited = set(["(0, 0)", "(7, 7)", "(4, 4)"])
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
    rootdir = '/home/DAVIDSON/doflocco/Graph_Pebbling/results/lemke square/ls-bounds'
    size = 16
    length = 16
    for r in nodes:
        if r not in visited: 
            strr = str(r).replace(" ", "")
            dir = "/v-" + strr
            dest = rootdir+dir
            # if not os.path.isfile(dest):
            #     r_rev = "("+ r[4]+ "," + r[1] +")"
            #     dest = rootdir + "/v-" + r_rev

            r1, r2 = int(r[1]), int(r[4])
            
            ls = PebblingGraph(lemkeSquare.edges, (r1, r2))
            optimizer = Optimizer(solver = "Gurobi", 
                                threads = None, 
                                logFile = "ls_sym_v-" + str(r),
                                lpFile = "ls_sym_v-" + str(r),
                                paramFile = "ls_sym_params" + "-v" + str(r),
                                objGap=False,
                                timeLimit=10800, 
                                cert=True)
            
            res = optimizer.generateTreeStrategiesSym(ls, size, length, None)
            #bounds[r] = res[1]
            strategies = res[0]
            results = res[2] 
            certDest = dest + "/certificates"
            #os.makedirs(os.path.dirname(certDest), exist_ok=True)
            edgeDest = dest + "/edges"
            #os.makedirs(os.path.dirname(edgeDest), exist_ok=True)

            for t in strategies.keys():
                # strategies[t].saveCertificate(certDest + "/ls_sym_cert"+ str(t)+ "-v" + str(r) +".csv")
                # strategies[t].saveEdges(edgeDest + "/ls_sym_edges" + str(t)+ "-v" + str(r) + ".csv")
                strategies[t].saveCertificate("ls_sym_cert"+ str(t)+ "-v" + str(r) +".csv")
                strategies[t].saveEdges("ls_sym_edges" + str(t)+ "-v" + str(r) + ".csv")
                #strategies[t].visualizeStrategy("ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            stratDest = dest 
            # os.makedirs(os.path.dirname(stratDest), exist_ok=True)
            # visualize strategy 
            #for t in strategies.keys():
            #     strategies[t].visualizeStrategy(stratDest + "/ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            #    strategies[t].visualizeStrategy("ls_sym_strategy" + str(t) + "-v" + str(r) + ".png")
            allResults = allResults.append(results, ignore_index=True)
            allResults.to_csv("results.csv")
        visited.add(r)
        r_rev = "("+ r[4]+ ", " + r[1] +")"
        visited.add(r_rev)
    
    pd.DataFrame.from_dict(bounds, orient="index").to_csv("lemke_square_bounds2.csv")


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

def readNetwork():
    
    dir = os.getcwd()
    ready = False

    while not ready: 
        print("SNDlib Network Library: ")
        print("[abilene, atlanta, cost226, eu, france, geant, germany, india, ny, ta2, us-ca, zib54]")
        name = input('Enter Network to Analyze: ')
   
        if name == 'abilene': 
            path = dir + '/networks/abilene/abilene--D-B-M-N-C-A-N-N.gml'
            ready = True
        elif name == 'alt' or name == 'atlanta': 
            path = dir + "/networks/atlanta/atlanta--D-B-M-N-C-A-N-N.gml"
            ready = True
        elif name == 'cost226':
            path = dir + "/networks/cost226/cost266--D-B-E-N-C-A-N-N.gml"
            ready = True
        elif name == 'eu' or name == 'europe':
            path = dir + '/networks/eruope/nobel-eu--D-B-E-N-C-A-N-N.gml'
            ready = True
        elif name == 'france': 
            path = dir + '/networks/france/france--D-B-L-N-C-A-N-N.gml'
            ready = True
        elif name == 'geant': 
            path = dir + '/networks/geant/geant--D-B-M-N-C-A-N-N.gml'
            ready = True
        elif name == 'germany': 
            path = dir + '/networks/germany/germany50--D-B-L-N-C-A-N-N.gml'
            ready = True
        elif name == 'india': 
            path = dir + '/networks/india/india35--D-B-M-N-C-A-N-N.gml'
            ready = True
        elif name == 'ny':
            path = dir + '/networks/new york /newyork--D-B-E-N-C-A-N-N.gml'
            ready = True
        elif name == 'ta2':
            path = dir + '/networks/ta2/ta2--D-B-E-N-C-A-N-N.gml'
            ready = True
        elif name == 'us-ca':
            path = dir + '/networks/us-ca/janos-us-ca--D-D-L-N-C-A-N-N.gml'
            ready = True
        elif name == 'zib54':
            path = dir + '/networks/zib54/zib54--U-U-E-N-I-A-N-N.gml'
            ready = True
        elif name == 'quit':
            exit(1)
        else: 
            print("Invalid Network Name, try again. Enter 'quit' to exit.")
    name = name + "_network"
    return loadNetwork(path), name
    

def loadNetwork(path):
    G = nx.read_gml(path)
    edges = list(set(G.edges))
    nodes = list(set(G.nodes))
    r = nodes[0]
    pebblingNetwork = PebblingGraph(edges, r)

    return pebblingNetwork

def getNetworkStats(pebblingNx):
    numNodes = len(pebblingNx.nodes)
    numEdges = len(pebblingNx.edges)

    nodeDegrees = list(pebblingNx.graph.degree(pebblingNx.nodes))
    maxDeg = max([d for v, d in nodeDegrees])
    minDeg = min([d for v, d in nodeDegrees])
    avgDeg = sum([d for v, d in nodeDegrees])/numNodes
    density = nx.density(pebblingNx.graph)
    diameter = nx.diameter(pebblingNx.graph)

    stats = {"|V|": numNodes, 
             "|E|": numEdges, 
             "Max degree": maxDeg,
             "Min degree": minDeg, 
             "Avg degree": avgDeg, 
             "Density": density, 
             "Diam": diameter}
    return stats

def applyPebbling():

    networkLib = ["abilene", "atlanta", "cost226", "eu", "france", "geant", "germany", "india", "ny", "ta2", "us-ca", "zib54"]
    #networkLib = ["us-ca", "germany", "cost226"]
    #headers = ["Tree Type", "Threads", "Solver", "|V|", "|E|", "|T|", "M"]
     
    allResults = pd.DataFrame()
    for network in networkLib:
        path = readNetwork(network)
        pebblingNetwork = loadNetwork(path)
        #print([item for item, count in collections.Counter(pebblingNetwork.arcs).items() if count > 1])
        nxStats = getNetworkStats(pebblingNetwork)

        
        size = 6
        length = 20 
        numNodes = None
        bounds = []
        for r in pebblingNetwork.nodes:
            results = {}
            optimizer = Optimizer(solver = "Gurobi", 
                                threads = None, 
                                logFile = network + "_v-" + str(r),
                                lpFile = network + "_v-" + str(r),
                                paramFile = network + "_v-" + str(r),
                                objGap=False,
                                timeLimit=1200, 
                                cert=False)
            pebblingNetwork.root = r
            res = optimizer.generateTreeStrategies(pebblingNetwork, size, length, numNodes)

            strategies = res[0]
            bound = res[1]
            results = res[2] 

            results["Network"] = network 
            results["Max degree"] = nxStats['Max degree']
            results["Min degree"] = nxStats['Min degree']
            results['Avg degree'] = nxStats['Avg degree']
            results["Density"] = nxStats['Density']

            bounds.append(bounds)

            for t in strategies.keys():
                # strategies[t].saveCertificate(network + "-v" + str(r) + "_cert" + str(t) +".csv")
                strategies[t].saveEdges(network + "_v-" + str(r) + "_edges-" + str(t) +".csv")
                strategies[t].visualizeStrategy(network + "_v-" + str(r) + "_strategy-" + str(t) + ".png")
            allResults = allResults.append(results, ignore_index=True)
            allResults.to_csv("network_pebbling_results.csv")
        # pebblingBound = max(bounds)


def loadGraph(name): 
    ready = False
    while not ready: 
        if name == "bruhat" or name == "Bruhat": 
            bruhat = [(0,1), (0,2), (0,4), (1,2), (1,7), (2,3), (2,20), (3,23), (4, 5), (4, 8), (5,6), (5,9), (6,7), (6,10), (7,11),
                (8,9), (8,16), (9,12), (10,11), (10,13), (11, 19), (12, 13), (12,14), (13,15), (14, 15), (14,17), (15, 18),
                (16, 17), (16,20), (17, 21), (18,19), (18,22), (19,23), (20, 21), (21,22), (22,23)]
            r = 0
            pebblingGraph = PebblingGraph(bruhat, r)
            ready = True
        elif name == "lemke" or name == "Lemke": 
            lemke_edges = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
            r = 0
            pebblingGraph = PebblingGraph(lemke_edges, r)
            ready = True
        elif name == "Petersen" or name == "petersen": 
            petersen = [(0,1), (1,2), (2,3), (3,4), (4,0), (0,6), (1,7), (2,8), (3,9), (4,5), (5,7), (5,8), (6,9), (6,8), (7,9)]
            r = 0
            pebblingGraph = PebblingGraph(petersen, r)
            ready = True
        elif name == "n-cube": 
            n = input("Enter value for n: ")
            nCube = nx.hypercube_graph(int(n))
            nCubeEdges = [e for e in nCube.edges]
            pebblingGraph = PebblingGraph(nCubeEdges, 0)
            ready = True
        elif name == "lemke square" or name == "Lemke Square" or name == "Lemke square": 
            lemke_edges = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
            r = 0
            lemke = PebblingGraph(lemke_edges, r)
            lemkeSquare = nx.cartesian_product(lemke.graph, lemke.graph)
            lemkeSquare = list(lemkeSquare.edges())
            pebblingGraph = PebblingGraph(lemkeSquare, r)
            ready = True

        elif name == "apply" or name == "Apply": 
            ready = True
            return readNetwork()
        elif name == 'quit':
            quit()
        else: 
            name = input("Invalid graph selection, please try again. Options are: \n [lemke, petersen, n-cube, bruhat, lemke square, network]\n>>> ")
            # need to make this loop back to initial prompt
    
    return pebblingGraph, name

def initializeOptimizer(fileRoot): 

    print("Default optimizer parameters: \n",
        "Solver : Gurobi \n", 
        "Threads : full \n", 
        "Root Filename :" + fileRoot +"\n", 
        "Objective Gap : None \n", 
        "Time Limit : None \n", 
        "Save Certificates? : True")
    change = input("Would you like to change optimizer paramters? y or n: " )
    if change == "n" or change == 'N':
        
        optimizer = Optimizer(solver="Gurobi", 
                              threads=False, 
                              logFile=fileRoot, 
                              lpFile=fileRoot,
                              paramFile=f"{fileRoot}_params",
                              timeLimit=False, 
                              objGap=False,
                              cert=True)
    else: 

        solver = input('Enter solver (Gurobi or CPLEX): ')
        threads = input('Enter number of threads (int): ')
        fileRoot = input("Enter root filename (str): ")
        timeLimit = input("Enter time limit (ms): ")
        objGap = input("Enter objective gap (float btwn 0 and 1): ")
        cert = input("Would you like to save certificate? (True or False): ")

        optimizer = Optimizer(solver=solver, 
                              threads=int(threads), 
                              logFile=fileRoot, 
                              lpFile=fileRoot,
                              paramFile=f"{fileRoot}_params",
                              timeLimit=int(timeLimit), 
                              objGap=float(objGap),
                              cert=cert)

    return optimizer

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
    r = (3,3)
    #print(countElementaryCirtuits(64))
    lemkeSquare = PebblingGraph(lemkeSquare, r)
    # treeStrategyOpt(lemkeSquare, 6, 15)
    # trimStrategiesSym(lemkeSquare)
    #print(lemkeSquare.nodes)
    #print(lemkeSquare.edges)
    # G = nx.read_gml('/home/DAVIDSON/doflocco/Graph_Pebbling/networks/us-ca/janos-us-ca--D-D-L-N-C-A-N-N.gml')
    # edges = [e for e in G.edges]
    # nodes = [v for v in G.nodes]
    # r = nodes[0]
    # pebblingNetwork = PebblingGraph(edges, r)
    # print(len(pebblingNetwork.nodes))
    # print(len(set(pebblingNetwork.nodes)))
    
    # applyPebbling()
    #treeStrategyOptCrank(lemkeSquare, 6, 12, None)
    boundLemkeSquare(lemkeSquare)
    #treeStrategyOptSymCrank(lemkeSquare, 20, 16, None)
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
    # networkLib = ["abilene", "atlanta", "cost226", "eu", "france", "geant", "germany", "india", "ny"]
    # results = pd.DataFrame()
    # for network in networkLib:
    #     print(network) 
    #     path = readNetwork(network)
    #     pebblingNetwork = loadNetwork(path)
    #     for r in pebblingNetwork.nodes: 
    #         file = network + "_v-" + r + ".csv"
    #         df = pd.read_csv(file)
    #         df['Network'] = [network]
    #         results = results.append(df, ignore_index=True)

        
    # results.to_csv("network_results.csv")

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
def setRoot(pebblingGraph): 
    ready = False
    while not ready: 
        print("Specify root from the node list:")
        print(pebblingGraph.nodes)
        print("To compute bound for all roots enter \"all\".")
        r = input(">>>> ")

        if r in pebblingGraph.nodes: 
            pebblingGraph.root = str(r)
            ready = True
        elif r == 'quit':
            quit()
        else: 
            print("Specified root vertex not in graph. Please enter a valid root from node list:")
            print(pebblingGraph.nodes)
            print("Or enter \"quit\" to exit program.")
    
    return pebblingGraph
def interface():

    print("Welcome to the Graph Pebbling Interactive MILP Platform!\n ")
    print("To begin, select your pebbling graph from the following preloaded library: ")
    graphName = input("[Lemke, Bruhat, n-cube, Petersen, Lemke square, apply (SNDlib network)] \n "
         ">>>> ")
    
    graphRes = loadGraph(graphName)
    pebblingGraph = graphRes[0]
    graphName = graphRes[1]
    print()
    print("Successfully Loaded %s graph." % graphName)
    print("---------------------------")
    print("Pebbling Graph Statistics: ")
    pebblingGraph.printStats() 
    print("---------------------------")
    print()

    pebblingGraph = setRoot(pebblingGraph)
    print("Specified root vertex for %s graph set to r = %s" % (graphName, pebblingGraph.root))
    print()
    
    print("Pebbling Graph ready for tree strategy generation.")
    print()

   
    fileRoot = f"{graphName}_r={pebblingGraph.root}"
    opt = initializeOptimizer(fileRoot)
    print()
    print("Optimization object successfully initialized with paramters: ")
    print("---------------------------")
    print(f"Solver : {opt.solver}\n",
          f"Threads : {opt.threads}\n", 
          f"Root Filename : {opt.logFile}\n",
          f"Time Limit : {opt.timeLimit}\n", 
          f"Objective Gap : {opt.objGap}\n", 
          f"Save Certificates? : {opt.cert}") 
    print("---------------------------")
    print()

    sym = input("Tree strategy generation (TS) or Symmetric Tree Strategy Generation (STS)? ")

    if sym == "ts" or sym == "TS": 
        size = input("Enter number of tree strategies to generate: ")
        length = input("Enter maximum tree length: ")
        print(f"Ready to generate {size} tree strategies with maximum tree length {length}")
        
    else: 
        size = input("Enter number of tree strategies to generate: ")
        length = input("Enter maximum tree length: ")
        print(f"Ready to generate {int(size)//2} symmetric tree strategies ({size} total) with maximum tree length {length}")
    input("Press Enter to run solver...")   
    print("------------------------------------------------------------------") 
    res = opt.generateTreeStrategies(pebblingGraph, int(size), int(length), None)
    
    strategies = res[0]
    if opt.cert: 
        for t in strategies.keys():
            # strategies[t].saveCertificate(f"{fileRoot}_cert{t}.csv")
            strategies[t].saveEdges(f"{fileRoot}_edges{t}.csv")

        #visualize strategy 
        for t in strategies.keys():
            strategies[t].visualizeStrategy(f"{fileRoot}_strategy{t}.png")

if __name__ == "__main__":
    interface()