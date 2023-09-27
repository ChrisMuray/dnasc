#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import json

T = 50
#Some constants
P_L = 2000                          # length of plasmid
TL_S = 203                          # length of transcribed portion of mSpinach
# TL_S = 68
# TL_G = 68
TL_G = 750
# PL_S = 40                           # length of p_Lac promoter
PL_S = 44                           # length of p_Lac promoter
PL_G = 44
# N_S = 150                           # length of intergenic region
N_S = 150

h0 = 10.5                           # number base pairs per turn in relaxed state
sigma0 = -1/h0                       # relaxed global supercoiling density
sigma_starS = sigma0 * P_L / PL_S  # optimal local supercoiling density
sigma_starG = sigma0 * P_L / PL_G

# Affect k_f and k_cat
zeta = 1 # elongation rate
beta = 1 # initiation rate

# Initialize system
s = {}
s['LacI']     = [500]
s['TetR']     = [500]
s['IPTG']     = [1000]
s['aTc']      = [1000]
s['aLacI']    = [0]
s['aTetR']    = [0]
s['pLacC']    = [0]
s['pTetC']    = [0]
s['p_Lac']    = [1000]
s['p_Tet']    = [1000]
s['R']        = [1000]
s['CC_S']     = [0]
s['EC_G']     = [0]
s['CC_G']     = [0]
s['EC_S']     = [0]
s['mS']       = [0]
s['MG']       = [0]
s['sigma_pS'] = [sigma0]
s['sigma_tS'] = [sigma0]
s['sigma_pG'] = [sigma0]
s['sigma_tG'] = [sigma0]

rho_L = 10
rho_T = 10
delta_p = 1
k_aL = 10
k_aT = 10
k_uaL = 1
k_uaT = 1
k_seqL = 10
k_seqT = 10
k_uL = 1
k_uT = 1
k_r = 1
k_open = 10
delta_m = 1
k_Mt = 1
k_Mg = 1
k_Msigma = 1
k_fS = zeta / (1 + (s['sigma_pS'][-1] - sigma_starS)**2 / k_Msigma)
k_fG = zeta / (1 + (s['sigma_pG'][-1] - sigma_starG)**2 / k_Msigma)
k_catS = beta / (1 + (s['sigma_tS'][-1] - sigma_starS)**2 / k_Msigma)
k_catG = beta / (1 + (s['sigma_tG'][-1] - sigma_starG)**2 / k_Msigma)

reactions = [
    [[],['LacI']],
    [['LacI'],[]],
    [['LacI','IPTG'],['aLacI']],
    [['aLacI'],['LacI','IPTG']],
    [['p_Lac','LacI'],['pLacC']],
    [['pLacC'],['p_Lac','LacI']],
    [['pLacC','IPTG'],['p_Lac','aLacI']],
    [['R','p_Lac'],['CC_S']],
    [['CC_S'],['R','p_Lac']],
    [['CC_S'],['EC_S']],
    [['EC_S'],['mS','R','p_Lac']],
    [['mS'],[]],
    [[],['TetR']],
    [['TetR'],[]],
    [['TetR','aTc'],['aTetR']],
    [['aTetR'],['TetR','aTc']],
    [['p_Tet','TetR'],['pTetC']],
    [['pTetC'],['p_Tet','TetR']],
    [['pTetC','aTc'],['p_Tet','aTetR']],
    [['R','p_Tet'],['CC_G']],
    [['CC_G'],['R','p_Tet']],
    [['CC_G'],['EC_G']],
    [['EC_G'],['MG','R','p_Tet']],
    [['MG'],[]]
]

def copy_system():
    for key in s.keys():
        s[key].append(s[key][-1])

def react(reactants, products):
    for r in reactants:
        s[r][-1] -= 1
    for p in products:
        s[p][-1] += 1

# Topoisomerase and gyrase dynamics
T0 = 10
topo_tau = 0.25
G0 = 10
gamma = 0.3
def m(sigma):
    d = sigma - sigma0
    if d > 0:
        return -G0 * gamma * d**2 / (k_Mg * sigma0 + d**2)
    else :
        return T0 * topo_tau * d**2 / (k_Mt * sigma0 + d**2)

t = [0]

reaction_indices = []
while t[-1] < T:

    # Update dynamic rates
    k_fS = zeta / (1 + (s['sigma_pS'][-1] - sigma_starS)**2 / k_Msigma)
    k_fG = zeta / (1 + (s['sigma_pG'][-1] - sigma_starG)**2 / k_Msigma)
    k_catS = beta / (1 + (s['sigma_tS'][-1] - sigma_starS)**2 / k_Msigma)
    k_catG = beta / (1 + (s['sigma_tG'][-1] - sigma_starG)**2 / k_Msigma)

    #Compute propensities
    props = np.asarray([
        rho_L,
        delta_p * s['LacI'][-1],
        k_aL * s['LacI'][-1] * s['IPTG'][-1],
        k_uaL * s['aLacI'][-1],
        k_seqL * s['p_Lac'][-1] * s['LacI'][-1],
        k_uL * s['pLacC'][-1],
        k_aL * s['pLacC'][-1] * s['IPTG'][-1],
        k_fS * s['R'][-1] * s['p_Lac'][-1],
        k_r * s['CC_S'][-1],
        k_open * s['CC_S'][-1],
        k_catS * s['EC_S'][-1],
        delta_m * s['mS'][-1],
        rho_T,
        delta_p * s['TetR'][-1],
        k_aT * s['TetR'][-1] * s['aTc'][-1],
        k_uaT * s['aTetR'][-1],
        k_seqT * s['p_Tet'][-1] * s['TetR'][-1],
        k_uT * s['pTetC'][-1],
        k_aT * s['pTetC'][-1] * s['aTc'][-1],
        k_fG * s['R'][-1] * s['p_Tet'][-1],
        k_r * s['CC_G'][-1],
        k_open * s['CC_G'][-1],
        k_catG * s['EC_G'][-1],
        delta_m * s['MG'][-1]
    ])

    prop_sum = sum(props)
    if prop_sum == 0:
        break

    tau = np.random.exponential(1/prop_sum)
    t.append(t[-1] + tau)

    r = np.random.rand() * prop_sum

    copy_system()

    for i in range(len(props)-1, -1, -1):
        if r > sum(props[:i]):
            reaction_indices.append(i)
            react(reactions[i][0], reactions[i][1])
            break 

    delta_kink = h0 / (s['sigma_pS'][-1] + s['sigma_pG'][-1])

    n_fS = max(PL_S + N_S/2 - delta_kink, 1)
    n_fG = max(PL_S + N_S/2 + delta_kink, 1)

    s['sigma_tS'][-1] +=  (PL_S / (2 * h0 * TL_S)) * (s['EC_S'][-1] - s['EC_S'][-2] - s['CC_S'][-1] + s['CC_S'][-2]) + m(s['sigma_tS'][-1])
    s['sigma_tG'][-1] +=  (PL_G / (2 * h0 * TL_G)) * (s['EC_G'][-1] - s['EC_G'][-2] - s['CC_G'][-1] + s['CC_G'][-2]) + m(s['sigma_tG'][-1])
    s['sigma_pS'][-1] += -(PL_S / (2 * h0 * n_fS)) * (s['EC_S'][-1] - s['EC_S'][-2] - s['CC_S'][-1] + s['CC_S'][-2]) + m(s['sigma_pS'][-1])
    s['sigma_pG'][-1] +=  (PL_G / (2 * h0 * n_fG)) * (s['EC_G'][-1] - s['EC_G'][-2] - s['CC_G'][-1] + s['CC_G'][-2]) + m(s['sigma_pS'][-1])

# reaction_didnt_happen = False
# for i in range(len(reactions)):
#     if (i not in reaction_indices):
#         if not reaction_didnt_happen:
#             print("Reactions that did not occur:")
#             reaction_didnt_happen = True
#         print(reactions[i])

# for k, v in s.items():
#     plt.plot(t, v, label=k)

# plt.subplot(1,2,1)
# plt.plot(t, s['sigma_tS'], label='sigma_tS')
# plt.plot(t, s['sigma_tG'], label='sigma_tG')
# plt.plot(t, s['sigma_pS'], label='sigma_pS')
# plt.plot(t, s['sigma_tG'], label='sigma_pG')
# plt.xlabel('Time (s)', fontsize=18)
# plt.ylabel('Supercoiling Density', fontsize=18)
# plt.title('Supercoiling')
# plt.legend()

# plt.subplot(1,2,2)
# plt.plot(t, s['mS'], label='mS')
# plt.plot(t, s['MG'], label='MG')
# plt.xlabel('Time (s)', fontsize=18)
# plt.ylabel('#', fontsize=18)
# plt.title('Protein Production')

# plt.legend()
# plt.savefig('plot.png')
