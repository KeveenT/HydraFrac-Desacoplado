#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:39:10 2020

@author: keveent
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

u_p = []
Pf = []

u_u = []
P = []

x_u = []
x_p = []
K = []

with open('poiseuille_u.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        u_p.append(row)

with open('poiseuille_p.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        Pf.append(row)

with open('ns_u.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        u_u.append(row)
    
with open('ns_p.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        P.append(row)

with open('x_p.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        x_p.append(row)

with open('x_u.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        x_u.append(row)

with open('K.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        K.append(row)

# x_p = np.array(x_p, dtype=float)
# u_p = np.array(u_p, dtype=float)
# Pf = np.array(Pf, dtype=float)
# x_u = np.array(x_u, dtype=float)
# u_u = np.array(u_u, dtype=float)
# L = np.array(L, dtype=float)
# P = np.array(P, dtype=float)

u_p = np.array(u_p, dtype=float)
Pf = np.array(Pf, dtype=float)
u_u = np.array(u_u, dtype=float)
P = np.array(P, dtype=float)

x_p = np.array(x_p, dtype=float)
x_u = np.array(x_u, dtype=float)
K = np.array(K, dtype=float)

# def get_normP():
#     Pf_norm = np.zeros((len(Pf), len(K)))
#     for j in range(0, len(K)):
#         for i in range(0, len(Pf)):
#             Pf_norm[i][j] = Pf[i][j] / Pf[0][j]
#     P_norm = np.zeros((len(P), len(K)))
#     for j in range(0, len(K)):
#         for i in range(0, len(P)):
#             P_norm[i][j] = P[i][j] / P[0][j]
#     return Pf_norm, P_norm

def get_normu():
    u_p_norm = np.zeros((len(u_p), len(K)))
    for j in range(0, len(K)):
        for i in range(0, len(Pf)):
            u_p_norm[i][j] = u_p[i][j] / u_p[0][j]
    u_u_norm = np.zeros((len(u_u), len(K)))
    for j in range(0, len(K)):
        for i in range(0, len(P)):
            u_u_norm[i][j] = u_u[i][j] / u_u[0][j]
    return u_p_norm, u_u_norm

# Pf_norm, P_norm = get_normP()
u_p_norm, u_u_norm = get_normu()

plt.ticklabel_format(axis="y", useOffset=False, style="sci")
# plt.ylim(-0.06, 1.6)
# plt.xlim(-0.3, 10.3)
# plt.text(8.1, 1.5, "K = 1e-15", size=10) #triangle
# plt.plot([5.6, 8.05], [0.98, 1.52], '-', color='k', linewidth=0.8)
# plt.text(8.1, 1.4, "K = 1e-14", size=10)
# plt.plot([5.4, 8.05], [0.85, 1.42], '-', color='k', linewidth=0.8)
# plt.text(8.1, 1.3, "K = 1e-13", size=10)
# plt.plot([4.4, 8.05], [0.517, 1.32], '-', color='k', linewidth=0.8)
# plt.text(8.1, 1.2, "K = 1e-12", size=10)
# plt.plot([3.3, 8.05], [0.18, 1.22], '-', color='k', linewidth=0.8)
# plt.text(8.1, 1.1, "K = 1e-11", size=10)
# plt.plot([3, 8.05], [0.01, 1.12], '-', color='k', linewidth=0.8)

plt.ylim(-60, 1600)
plt.xlim(-0.3, 11)
# plt.text(8.8, 1390, "K = 1e-16", size=10) #triangle
# plt.plot([7.58, 8.750], [990, 1410], '-', color='k', linewidth=0.8)
plt.text(8.8, 1270, "K = 1e-15", size=10)
plt.plot([7.76, 8.75], [940, 1300], '-', color='k', linewidth=0.8)
plt.text(8.8, 1150, "K = 1e-14", size=10)
plt.plot([7.3, 8.750], [680, 1200], '-', color='k', linewidth=0.8)
plt.text(5.8, 1390, "K = 1e-13", size=10)
plt.plot([3.8, 5.750], [570, 1430], '-', color='k', linewidth=0.8)
plt.text(5.8, 1270, "K = 1e-12", size=10)
plt.plot([3.22, 5.750], [170, 1300], '-', color='k', linewidth=0.8)
plt.text(5.8, 1150, "K = 1e-11", size=10)
plt.plot([3.2, 5.750], [10, 1170], '-', color='k', linewidth=0.8)

# plt.ylim(-0.06, 1.6)
# plt.xlim(-0.3, 11)
# plt.text(8.8, 1.39, "K = 1e-16", size=10) #square
# plt.plot([7.58, 8.75], [0.99, 1.41], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.27, "K = 1e-15", size=10)
# plt.plot([7.9, 8.75], [0.99, 1.29], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.15, "K = 1e-14", size=10)
# plt.plot([8, 8.75], [0.91, 1.17], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.39, "K = 1e-13", size=10)
# plt.plot([3.85, 5.75], [.62, 1.41], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.27, "K = 1e-12", size=10)
# plt.plot([3.23, 5.75], [0.22, 1.29], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.15, "K = 1e-11", size=10)
# plt.plot([3, 5.75], [0.01, 1.17], '-', color='k', linewidth=0.8)

# plt.ylim(-0.06, 1600)
# plt.xlim(-0.3, 11)
# plt.text(8.8, 1390, "K = 1e-16", size=10) #square
# plt.plot([7.58, 8.75], [990, 1410], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1270, "K = 1e-15", size=10)
# plt.plot([7.9, 8.75], [990, 1290], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1150, "K = 1e-14", size=10)
# plt.plot([8, 8.75], [910, 1170], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1390, "K = 1e-13", size=10)
# plt.plot([3.85, 5.75], [620, 1410], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1270, "K = 1e-12", size=10)
# plt.plot([3.23, 5.75], [220, 1290], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1150, "K = 1e-11", size=10)
# plt.plot([3, 5.75], [10, 1170], '-', color='k', linewidth=0.8)

# plt.ylim(-0.06, 1600)
# plt.xlim(-0.3, 11)
# # plt.text(8.8, 1390, "K = 1e-16", size=10) #square
# # plt.plot([7.58, 8.75], [990, 1410], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1270, "K = 1e-15", size=10)
# plt.plot([7.9, 8.75], [990, 1290], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1150, "K = 1e-14", size=10)
# plt.plot([7.7, 8.75], [820, 1170], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1390, "K = 1e-13", size=10)
# plt.plot([3.85, 5.75], [620, 1410], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1270, "K = 1e-12", size=10)
# plt.plot([3.23, 5.75], [220, 1290], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1150, "K = 1e-11", size=10)
# plt.plot([3, 5.75], [10, 1170], '-', color='k', linewidth=0.8)



plt.plot(x_p, [row[0] for row in Pf], marker='.', markerfacecolor='k', markeredgecolor='k', markevery=int(len(x_p)/25), markersize=6, linestyle='None', label='Poiseuille')
plt.plot(x_p, [row[0] for row in P], color='k', linewidth=1.2, label='Navier-Stokes')
for i in range(1, len(K)):
    plt.plot(x_p, [row[i] for row in Pf], marker='.', markerfacecolor='k', markeredgecolor='k', markevery=int(len(x_p)/25), markersize=6, linestyle='None',)
    plt.plot(x_p, [row[i] for row in P], color='k', linewidth=1.2)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Pressão [Pa]")
# plt.title("Comparação de Pressões")
plt.legend(loc=2)
plt.grid()
plt.savefig('press_triangle_40_pp.pdf', dpi=1200)
plt.figure()


plt.ylim(-0.06, 1.6)
plt.xlim(-0.3, 11)
# plt.text(8.8, 1.39, "K = 1e-16", size=10) #triangle
# plt.plot([7.55, 8.75], [0.95, 1.41], '-', color='k', linewidth=0.8)
plt.text(8.8, 1.27, "K = 1e-15", size=10)
plt.plot([7.65, 8.75], [0.85, 1.29], '-', color='k', linewidth=0.8)
plt.text(8.8, 1.15, "K = 1e-14", size=10)
plt.plot([7.1, 8.75], [0.54, 1.17], '-', color='k', linewidth=0.8)
plt.text(5.8, 1.39, "K = 1e-13", size=10)
plt.plot([3.6, 5.75], [0.5, 1.41], '-', color='k', linewidth=0.8)
plt.text(5.8, 1.27, "K = 1e-12", size=10)
plt.plot([3.1, 5.75], [0.16, 1.29], '-', color='k', linewidth=0.8)
plt.text(5.8, 1.15, "K = 1e-11", size=10)
plt.plot([3, 5.75], [0.01, 1.17], '-', color='k', linewidth=0.8)

# plt.ylim(-0.06, 1.6)
# plt.xlim(-0.3, 11)
# plt.text(8.8, 1.39, "K = 1e-16", size=10) #square
# plt.plot([6.4, 8.75], [0.52, 1.41], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.27, "K = 1e-15", size=10)
# plt.plot([6.5, 8.75], [0.45, 1.29], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.15, "K = 1e-14", size=10)
# plt.plot([6.4, 8.75], [0.28, 1.17], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.39, "K = 1e-13", size=10)
# plt.plot([3.57, 5.75], [0.52, 1.41], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.27, "K = 1e-12", size=10)
# plt.plot([3.2, 5.75], [0.22, 1.29], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.15, "K = 1e-11", size=10)
# plt.plot([3, 5.75], [0.01, 1.17], '-', color='k', linewidth=0.8)

# plt.ylim(-0.06, 1.6)
# plt.xlim(-0.3, 11)
# plt.text(8.8, 1.39, "K = 1e-16", size=10) #square
# plt.plot([6, 8.75], [0.4, 1.41], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.27, "K = 1e-15", size=10)
# plt.plot([6.3, 8.75], [0.37, 1.29], '-', color='k', linewidth=0.8)
# plt.text(8.8, 1.15, "K = 1e-14", size=10)
# plt.plot([6.5, 8.75], [0.33, 1.17], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.39, "K = 1e-13", size=10)
# plt.plot([3.57, 5.75], [0.52, 1.41], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.27, "K = 1e-12", size=10)
# plt.plot([3.2, 5.75], [0.22, 1.29], '-', color='k', linewidth=0.8)
# plt.text(5.8, 1.15, "K = 1e-11", size=10)
# plt.plot([3, 5.75], [0.01, 1.17], '-', color='k', linewidth=0.8)

plt.plot(x_u, [row[0] for row in u_p_norm], marker='.', markerfacecolor='k', markeredgecolor='k', markevery=int(len(x_p)/25), markersize=6, linestyle='None', label='Poiseuille')
plt.plot(x_u, [row[0] for row in u_u_norm], color='k', linewidth=1.2, label='Navier-Stokes')
for i in range(1, len(K)):
    plt.plot(x_u, [row[i] for row in u_p_norm], marker='.', markerfacecolor='k', markeredgecolor='k', markevery=int(len(x_p)/25), markersize=6, linestyle='None',)
    plt.plot(x_u, [row[i] for row in u_u_norm], color='k', linewidth=1.2)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Velocidade Normalizada")
# plt.title("Comparação de Velocidades")
plt.legend(loc=2)
plt.grid()
plt.savefig('vel_triangle_40_pp.pdf', dpi=1200)
plt.figure()


plt.show()