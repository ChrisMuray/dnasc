#!/usr/bin/python3

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# print(sys.argv)

if len(sys.argv) != 3:
    print("Usage:", sys.argv[0], "[filename] [feature name]")
    exit()

infile = open(sys.argv[1])
feature_name = sys.argv[2]

sims = []
for i in range(1000):
    sims.append(json.loads(infile.readline()))
infile.close

max_T = 0
for s in sims:
    max_T = max(max_T, s['t'][-1])

num_timebins = 100
t = np.linspace(0, max_T, num_timebins + 1)
with open(feature_name + "_binned.json", "w") as outfile:
    json.dump(t.tolist(), outfile)
    outfile.write("\n")

features = []
for i in range(num_timebins):
    current_feature = []
    sim_num = 0
    for s in sims:
        # if sim_num >= 10:
        #     break
        for j in range(len(s['t'])):
            if s['t'][j] >= t[i] and s['t'][j] < t[i+1]:
                current_feature.append(s[feature_name][j])
        # sim_num += 1
    features.append(current_feature)

with open(feature_name + "_binned.json", "a") as outfile:
    json.dump(features, outfile)