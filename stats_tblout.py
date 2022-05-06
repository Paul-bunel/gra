import sys
import os
from statistics import *
import json

def stats(dbpath):
    total_results = []
    dir_res = {}
    for file in os.listdir(dbpath):
        filename = os.path.join(dbpath, file)
        results = []
        with open(filename, 'r') as f:
            for line in f:
                if line[0] != '#':
                    results.append(line.split())
                    total_results.append(line.split())
        dir_res[file] = [
            mean([float(row[4]) for row in results]),
            min([float(row[4]) for row in results]),
            max([float(row[4]) for row in results])
        ]

    dir_res["general_mean"] = mean([float(row[4]) for row in total_results])
    return dir_res


def main():
    HMMstats = {}
    if len(sys.argv) == 1:
        databases = os.listdir("RESULTS")
    else:
        databases = sys.argv[1:]

    for db in databases:
        dbpath = os.path.join("RESULTS", db)
        HMMstats[db] = stats(dbpath)

    print(f"HMM_PROFILE\t{'DB SCANNED':<30}   {'MEAN':<17}{'MIN':<10}{'MAX':<10}")
    print("_____________________________________________________________________________________\n")

    for k in HMMstats.keys():
        print(k, ':')
        [print(f"\t{key :<30}: {value}") for key, value in HMMstats[k].items()]


    json_str = json.dumps(HMMstats, skipkeys = False, indent = 4)
    with open("stats.json", 'w') as file:
        file.write(json_str)
        
main()
