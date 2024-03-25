# Automating Weight Function Generation in Graph Pebbling
Mixed-Integer Linear Programming Approaches to Weight Function Generation in Graph Pebbling

This code base provides an inquiry into Mixed-Integer Linear Programming (MILP) approaches for Weight Function generation in Graph Pebbling. The work is a culimation of an Applied Mathematics Honors Thesis at Davidson College under the supervision of Dr. Carl Yerger and Dr. Jonad Pulaj, and supplements the formal written report. This repository is supplemental material to the following manuscrtipt: 

> Dominic Flocco, Jonad Pulaj, and Carl Yerger. "Automating weight function generation in graph pebbling". _Discrete Applied Mathematics_ 347 (2024): 155-174.

In addition to including the certificates that prove pebbling bounds on the Bruhat and Lemke square graph, this code base in constructed for further experiementation by future researchers. The novel result of the thesis was a verifiable proof that the pebbling number of the Lemke square is upper bounded by 96. However, the computational framework provided is capable of computing pebbling upper bounds on general graphs, including real-world networks, some of which are preloaded into the codebase. 

## File Structure 

Optimizer.py -> mixed-integer linear program implementation <br/>
PebblingGraph.py -> pebbling graph class with appropriate methods for functionality<br/>
TreeStrategy.py -> weight function classes <br/>
main.py -> compiler file used to run experiments and test function classes <br/>

## Experiementation
The repository is structured as follows: <br/>
### Graph Library 
The following graphs are preloaded as Pebbling Graph Objects into the code base: 
- Petersen 
- Lemke
- n-cube 
- 4th order Bruhat 
- Lemke square <br>

We also include a number of netowrks compiled from [SNDlib](http://sndlib.zib.de/home.action): 
- Abilene 
- Atlanta
- Cost226 
- Europe 
- France 
- Geant 
- Germany 
- India 
- New York 
- US-CA 
- ta2

To view pebbling graph statistics (number of nodes, edges, diameter, max and min degree, etc.) load any of the above graphs into the compiler. 
### Running the MILP 

To run the MILP solver, compile and run the main.py file. This will prompt an interactive interface in the console that allows the user to load a pebbling graph, set MILP parameters and generate tree strategies. We provide workflow for the Lemke graph below. 

To begin, choose which graph you would like to pebble from the list above: 
```
To begin, select your pebbling graph from the following preloaded library: 
[Lemke, Bruhat, n-cube, Petersen, Lemke square, apply (SNDlib network)]
 >>>> lemke
````
When the graph is successfully loaded, the graphs statistics will be ouputted into the console. 

```
Successfully Loaded lemke graph.
---------------------------
Pebbling Graph Statistics: 
|V|: 8
|E|: 13
Max degree: 5
Min degree: 2
Avg degree: 3
Density: 0
Diam: 3
---------------------------
```

Next, select a root vertex from the graph's node list. 
```
Specify root from the node list:
['4', '3', '2', '0', '5', '7', '6', '1']
To compute bound for all roots enter "all".
>>>> 3
Specified root vertex for lemke graph set to r = 3

Pebbling Graph ready for tree strategy generation.
```

Now that the pebbling graph object is initialized with the specified root vertex, the optimizer object must be initialized. The console will begin by outputting the default parameters and give the user the option to change the parameters. 
```
Default optimizer parameters: 
 Solver : Gurobi 
 Threads : full 
 Root Filename :lemke_r=3
 Objective Gap : None 
 Time Limit : None 
 Save Certificates? : True
Would you like to change optimizer paramters? y or n: n

Optimization object successfully initialized with paramters: 
---------------------------
Solver : Gurobi
 Threads : False
 Root Filename : "lemke_r=3
" Time Limit : False
 Objective Gap : False
 Save Certificates? : True
---------------------------
```

Lastly, after initializing the optimizer object, the user must specify whether to run the Vanilla Tree Strategy Optimization Problem (TS) or the Symmetric Tree Strategy Optimization Problem (STS). Note that STS may only be run on Cartesian Product Graphs. In addition, the console will prompt the user to specify how many strategies to generate and to set the maximum tree length.  
```
Tree strategy generation (TS) or Symmetric Tree Strategy Generation (STS)? ts
Enter number of tree strategies to generate: 5
Enter maximum tree length: 6
```
To run the program, confirm that the paramters are correct and press enter. This will run the desired optimization software to generate tree strategies. The results, visualizations and log files will be saved into the current directory, and a summary of results will be printed to the console after the MILP has finished running. 

```
Ready to generate 5 tree strategies with maximum tree length 6
Press Enter to run solver...
-------------------------------------------------------------------
Academic license - for non-commercial use only - expires 2022-09-22
Using license file /Users/dominicflocco/gurobi.lic
Changed value of parameter LogFile to lemke_r=3.log
   Prev:   Default: 
Changed value of parameter MIPFocus to 1
   Prev: 0  Min: 0  Max: 3  Default: 0
Warning: variables 0 and 42 have the same name "X[2,4]"
Warning: linear constraint 45 and linear constraint 46 have the same name "weight-constr-t4"
Warning: to let Gurobi read it back, use rlp format
Gurobi Optimizer version 9.1.2 build v9.1.2rc0 (mac64)
Thread count: 2 physical cores, 4 logical processors, using up to 4 threads
Optimize a model with 212 rows, 210 columns and 630 nonzeros
Model fingerprint: 0x3a0e032d
Variable types: 40 continuous, 170 integer (170 binary)
Coefficient statistics:
  Matrix range     [1e+00, 6e+01]
  Objective range  [2e-01, 2e-01]
  Bounds range     [1e+00, 3e+01]
  RHS range        [1e+00, 6e+01]
Found heuristic solution: objective 12.4000000
Presolve removed 90 rows and 110 columns
Presolve time: 0.00s
Presolved: 122 rows, 100 columns, 390 nonzeros
Variable types: 30 continuous, 70 integer (65 binary)

Root relaxation: objective 7.000000e+00, 38 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    7.00000    0    6   12.40000    7.00000  43.5%     -    0s
H    0     0                       9.0000000    7.00000  22.2%     -    0s
     0     0    7.00000    0    4    9.00000    7.00000  22.2%     -    0s
H    0     0                       8.2000000    7.00000  14.6%     -    0s
     0     0    7.00000    0    9    8.20000    7.00000  14.6%     -    0s
     0     0    7.00000    0    2    8.20000    7.00000  14.6%     -    0s
H    0     0                       7.5000000    7.00000  6.67%     -    0s
     0     0    7.00000    0    6    7.50000    7.00000  6.67%     -    0s
     0     0    7.00000    0    7    7.50000    7.00000  6.67%     -    0s
     0     0    7.00000    0    4    7.50000    7.00000  6.67%     -    0s
H    0     0                       7.4000000    7.00000  5.41%     -    0s
     0     0    7.00000    0    6    7.40000    7.00000  5.41%     -    0s
H    0     2                       7.2000000    7.00000  2.78%     -    0s
     0     2    7.00000    0    6    7.20000    7.00000  2.78%     -    0s
H   11     3                       7.0000000    7.00000  0.00%   9.8    0s

Cutting planes:
  Learned: 1
  Cover: 1
  Clique: 1
  MIR: 6
  Flow cover: 5
  Relax-and-lift: 1

Explored 14 nodes (391 simplex iterations) in 0.08 seconds
Thread count was 4 (of 4 available processors)

Solution count 7: 7 7.2 7.4 ... 12.4

Optimal solution found (tolerance 1.00e-04)
Best objective 7.000000000000e+00, best bound 7.000000000000e+00, gap 0.0000%
----------------------------
Tree Type: symmetric
Threads: False
Solver: Gurobi
|V|: 8
|E|: 13
|T|: 5
Num Nodes: unbounded
Max Len: 6
runtime: 0.07922601699829102
root: 3
Bound: 8
----------------------------
```

## Pebbling Certificates 
Included in this repository are the weight function certificates and edges for the tree strategies that prove an upper bound of 96 on the Lemke square. The certificates for each vertex can be found in the directory results -> lemke square -> ls-bounds. 


<!-- ### visualizations 
saved netowrkx png vizualization files that display tree strategies generated from optimiation problems. (Lemke square visualizations are work in progress) 
### lps 
lp files for optimization problems run using Gurobi solver<br/>
### logs 
linear programming logs that track progress of optimization over time <br/>
### lemke_trees 
csv files that store edge sets for lemke square strategies generated by linear programs <br/>
### lemke_certificates 
weight function certificates for Lemke square generated by linear programs <br/> -->

