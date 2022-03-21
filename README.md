# Graph Pebbling
Mixed-Integer Linear Programming Approaches to Weight Function Generation in Graph Pebbling

This code base provides an inquiry into Mixed-Integer Linear Programming (MILP) approaches for 
Weight Function generation in Graph Pebbling. The work is a culimation of one semester of work 
in pursuit of an Applied Mathematics Honors Thesis at Davidson College under the supervision of 
Dr. Carl Yerger and Dr. Jonad Pulaj, and supplements a formal writen final report for the Fall 2021
semester. 
## Repository Structure
The repository is structured as follows: <br/>
### py 
Optimizer.py -> mixed-integer linear program implementation <br/>
PebblingGraph.py -> pebbling graph class with appropriate methods for functionality<br/>
TreeStrategy.py -> weight function classes <br/>
main.py -> compiler file used to run experiments and test function classes <br/>
### visualizations 
saved netowrkx png vizualization files that display tree strategies generated from optimiation problems. (Lemke square visualizations are work in progress) 
### lps 
lp files for optimization problems run using Gurobi solver<br/>
### logs 
linear programming logs that track progress of optimization over time <br/>
### lemke_trees 
csv files that store edge sets for lemke square strategies generated by linear programs <br/>
### lemke_certificates 
weight function certificates for Lemke square generated by linear programs <br/>
### bruhat_exp (out dated) 
stores results of optimization experiment on max length and set size for 4th weak Bruhat<br/>
