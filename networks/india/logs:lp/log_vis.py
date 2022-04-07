import grblogtools as glt
import os
import pandas as pd
import glob

def main():
    #print(glt.get_dataframe(["/home/doflocco/optimization/honors_thesis/lemke-s30-l10.log"], timelines=True))
    cwd = os.getcwd()
    logfiles = [cwd+'/trim_ls-size='+i for i in range(4, 15)]
    print(logfiles)
    summary, timelines = glt.get_dataframe(["/home/doflocco/optimization/honors_thesis/lemke-s30-l10.log"], timelines=True)
    print(summary)
    glt.plot(summary, type="box")

    #glt.plot(timelines, y="Gap", color="Log", type="line")

if __name__ == "__main__":
    main()