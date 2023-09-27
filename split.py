import json
import sys

with open("sims.json") as infile:
    sims = json.load(infile)

for i in range(len(sims)):
    with open("sims_split.json","a") as outfile:
        outfile.dump(sims[i])
        outfile.write("\n")