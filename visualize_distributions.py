#!/usr/bin/python3

import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "[filename]")
    exit()

with open(sys.argv[1]) as infile:
    sims = json.load(infile)

max_T = 0
for s in sims:
    max_T = max(max_T, s['t'][-1])

print(max_T)

num_timebins = 100
t = np.linspace(0, max_T, num_timebins + 1)

mS = []
for i in range(num_timebins):
    current_mS = []
    for j in range(len(s['mS'])):
        if s['mS'][j] >= t[i] and s['mS'][j] < t[i+1]:
            current_mS.append(s['mS'][j])
    mS.append(np.asarray(current_mS))

os.mkdir('histograms')
for i in range(num_timebins):
    plt.hist(mS[i])
    plt.title('Overall mS distribution from t = ' + str(t[i]) + ' to ' + str(t[i+1]))
    plt.savefig('histograms/mS'+str(i)+'.png')