import os 
import shutil 

for file in os.listdir("/home/DAVIDSON/doflocco/Graph_Pebbling/py"):
    if "params" in file:
        new_dir = "/home/DAVIDSON/doflocco/Graph_Pebbling/results/lemke square/ls4-bounds/parameters/" + os.path.basename(file)
        os.rename(file, new_dir)