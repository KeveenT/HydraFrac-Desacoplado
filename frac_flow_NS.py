#Programa para a simulação de escoamento em fratura

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve
import csv
import datetime
begin_time = datetime.datetime.now()

#Propriedades
rho = 1000 #1/Pa
mu = 1e-3
K = 1e-12

#Geometria
L = 10 #Comprimento da Fratura
D0 = .001 #Abertura da entrada (em metros)
# D0 = 0.75
#Malha e Tempo
dt = 1
tf = 100000000000000000000
nt = int(tf / dt)

nx = 1000 #Nบmero de Volumes
nu = (nx + 1) #Nรบmero de volumes de velocidade
m = nu + nx
dx = L/nx
Ln = 0.1
Ls = 0.1
tol = 1e-6 #Tolerรขncia para diferenรงa entre iteraรงรตes de tempo para regime permanente

#Aberturas da Fratura
def get_D(D0):
    global D
    D = np.zeros(nx+1)
    for i in range(0, len(D)): #Abertura de Abertura Constante
        D[i] = D0 
    # D = np.linspace(D0, D0/100, num = int(len(D)), endpoint=True) #Abertura Triangular
    D = np.linspace(D0, D0/10, num = int(len(D)), endpoint=True)
    return D
    
def get_Dp():
    global Dp
    D = get_D(D0)
    Dp = np.zeros(nx)
    for i in range(0, len(Dp)):
        Dp[i] = ((D[i]+D[i+1])/2)
    return Dp

#Montagem de Matriz de Coeficientes
def get_Auup():
    D = get_D(D0)
    Dp = get_Dp()
    Auup = np.zeros(nu)
    Auup[0] = ((rho * dx * D[0]) / (dt)) + ((2 * mu * Dp[0]) / (dx)) + 2*((6 * dx * mu)/ (D[0]))
    if D[-1] == 0.0:
        Auup[-1] = 1
    else:
        Auup[-1] = ((rho * dx * D[-1]) / (dt)) + ((2 * mu * Dp[-1]) / (dx)) + 2*((6 * dx * mu)/ (D[-1]))
    for i in range(1, len(D)-1):
        Auup[i] = ((rho * dx * D[i]) / (dt)) + ((2 * mu * Dp[i]) / (dx)) + ((2 * mu * Dp[i-1]) / (dx)) + 2*((6 * dx * mu)/ (D[i])) 
    return Auup

def get_Auuw():
    Dp = get_Dp()
    Auuw = np.zeros(nu-1)
    for i in range(0, len(Auuw)):
        Auuw[i] = 2*(mu * Dp[i]) / (dx)
    return Auuw

def get_Auue():
    Dp = get_Dp()
    Auue = np.zeros(nu-1)
    for i in range(0, len(Auue)):
        Auue[i] = 2*(mu * Dp[i]) / (dx)
    return Auue

def get_Aupe():
    D = get_D(D0)
    Dp = get_Dp()
    Aupe = np.zeros(nx)
    for i in range(0, len(Aupe)):
        Aupe[i] = Dp[i]
    return Aupe

def get_Aupw():
    D = get_D(D0)
    Dp = get_Dp()
    Aupw = np.zeros(nx)
    for i in range(0, len(Aupw)):
        Aupw[i] = Dp[i]
    return Aupw

def get_Aupp():
    Aupp = np.zeros(nx)
    for i in range(0, len(Aupp)):
        Aupp[i] = ((rho * dx * K) / mu) * ((1/Ln) + (1/Ls))
    return Aupp

#Montagem de Matriz de Coeficientes
def build_A(): #Monta a matriz A
    A = np.zeros([m,m])
    Auup = get_Auup()
    Auuw = get_Auuw()
    Auue = get_Auue()
    Aupw = get_Aupw()
    Aupe = get_Aupe()
    Aupp = get_Aupp()
    for i in range(0, len(Auup)):
        A[i, i] = Auup[i]
    for i in range(0, len(Auue)):
        A[i, i+1] = -Auue[i]
    for i in range(0, len(Auuw)):
        A[i+1, i] = -Auuw[i]
    for i in range(1, nx): #Adiciona o gradiente de pressão à matriz de QM
        A[i, nx+i] = -Aupw[i]
        A[i, nx+i+1] = Aupe[i]
    for i in range(0, len(Aupp)):
        A[nu+i, nu+i] = Aupp[i]
    A[0, nu] = 2 * ((D[0] + Dp[0])/2) #Condições de Contorno de pressão prescrita
    # A[nx, nu+nx-1] = 2 * (-(D[-1] + Dp[-1])/2)
    # A[0, 0] = 1
    # for i in range(1, nx+nu):
    #     A[0, i] = 0
    for i in range(nx, nx+nu):
        A[nx, i-1] = 0
    A[nx, nx] = 1
    # for i in range(0, nx+nu):
    #     A[nu, i-1] = 0
    # A[nu, nu] = 1
    return A

def get_A(S_star): #Retorna a matriz A com os termos advectivos
    global A
    A = build_A()
    D = get_D(D0)
    Dp = get_Dp()
    ue = []
    uw = []
    for i in range(0, nu-1): #Velocidade da conservação da massa 
        ue.append((S_star[i]+S_star[i+1])/2)
    for i in range(1, nu):
        uw.append((S_star[i]+S_star[i-1])/2)
    for k in range(0, len(ue)): #Algoritmo para Upwind
        if ue[k] > 0 and uw[k] > 0:
            A[k, k] += rho * ue[k] *Dp[k] #Adição de fluxo de massa a submatriz da conservação de quantidade de movimento
            A[k+1, k] += -rho * uw[k] * Dp[k]
        else:
            A[k, k+1] += rho * ue[k] * Dp[k]
            A[k+1, k+1] += -rho * uw[k] * Dp[k]
    # A[nx, nx] += rho * ue[-1] *Dp[-1]
    for m in range(0, nx): #Adição dos fluxos de massa a submatriz da conservação da massa
        A[nu+m, m] += -rho * D[m]
    for m in range(0, nx):
        A[nu+m, m+1] += rho * D[m+1]
    A[0, nu] = 2 * ((D[0] + Dp[0])/2) #Condições de Contorno de pressão prescrita
    # A[0, 0] = 1
    # for i in range(1, nx+nu):
    #     A[0, i] = 0
    for i in range(nx, nx+nu):
        A[nx, i-1] = 0
    A[nx, nx] = 1
    # for i in range(0, nx+nu):
    #     A[nu, i-1] = 0
    # A[nu, nu] = 1
    return A

#Construção do Termo Fonte
def get_Pn():
    Pn = np.zeros((nx,1))
    for i in range(0, len(Pn)):
        Pn[i] = 0
    return Pn

def get_Ps():
    Ps = np.zeros((nx,1))
    for i in range(0, len(Ps)):
        Ps[i] = 0
    return Ps

def get_b(S_old):
    global P_entrada, P_saida
    D = get_D(D0)
    Dp = get_Dp()
    Pn = get_Pn()
    Ps = get_Ps()
    b = np.zeros((nu+nx,1))
    massa = 0.001*rho
    P_entrada = 5e6
    P_saida = 0
    for i in range(0, nu):
        b[i] += ((rho * dx * D[i] * S_old[i]) / (dt))
    for i in range(nu, nx+nu):
        b[i] += ((((rho * dx * K) / mu) * (1 / Ln)) * Pn[i-nu]) + ((((rho * dx * K) / mu) * (1 / Ls)) * Pn[i-nu])
    # b[0] += massa/(rho * D[0])
    # b[nx] = 0
    b[0] +=  2 * P_entrada * ((D[0] + Dp[0])/2)
    # b[-nx] += P_entrada
    # b[nu-1] += -2 * P_saida * ((D[-1] + Dp[-1])/2)
    return b

#Solver
def solve():
    D = get_D(D0)
    S_old = np.zeros((nu+nx,1))
    for i in range(nu, len(S_old)):
        S_old[i] = 0
    m_ponto = 1
    # for i in range(0, nu):
    #     S_old[i] = 0#m_ponto / (rho * D[i])
    S_star = np.copy(S_old)
    t = 0
    while t < tf:
        n = 0
        b = get_b(S_old)
        while n < 1:
            A = get_A(S_star)
            lu, piv = lu_factor(A)
            S_star = lu_solve((lu, piv), b)
            n = n + 1
        if (np.linalg.norm(S_star-S_old,np.inf) < tol):
            print(t)
            break
        else:
            S_old = S_star.copy()
            t = t + dt
    if np.allclose(np.dot(A, S_star), b) == True:
        print('Solução Encontrada')
    else:
        print('Não foi obtida uma solução')
    return S_star

def get_Slist():
    global K, Klist
    Slist = []
    Klist = []
    K = 1e-15
    while K <= 1e-15:
        print(K)
        Klist.append(K)
        S_star = solve()
        Slist.append(S_star)
        K = K*10
    return Slist

Slist = get_Slist()

def get_results_u(Slist):
    results_u = np.zeros((nu, len(Slist)))
    for h in range(0, len(Slist)):
        for i in range(0, nu):
            results_u[i][h] = Slist[h][i][0]
    return results_u

def get_results_p(Slist):
    results_p = np.zeros((nx, len(Slist)))
    for h in range(0, len(Slist)):
        for i in range(nu, nu+nx):
            results_p[i-nu][h] = Slist[h][i][0]
    return results_p

results_u = get_results_u(Slist)
results_p = get_results_p(Slist)
# S, b = solve()
# u = np.zeros((nu,1))
# P = np.zeros((nx,1))
# for i in range(0, nu):
#     u[i] = S[i]

# for j in range(0, len(u)):
#     u[j] = round(u[j][0], 6)
    
# for i in range(nu, nu+nx):
#     P[i-nu] = S[i]

# x = np.zeros((nu,1))
# for i in range(0, len(x)):
#     x[i] = x[i-1] + dx
    
# m_ponto = np.zeros(len(u))
# for i in range(0, len(u)):
#     m_ponto[i] = (rho * D0) * u[i]

# plt.plot(x, m_ponto)
# plt.xlabel("Comprimento [m]")
# plt.ylabel("m ponto [kg/s]")
# plt.title("m ponto - Navier-Stokes")
# plt.grid()
# plt.figure()

# L = np.zeros((nx,1))
# for i in range(0, len(L)):
#     L[i] = L[i-1] + dx 
# plt.plot(L, P)
# plt.xlabel("Comprimento [m]")
# plt.ylabel("Pressão [Pa]")
# plt.title("Pressão - Navier-Stokes")
# plt.grid()
# plt.figure()

# plt.ticklabel_format(axis="y", useOffset=False, style="plain")
# plt.plot(x, u)
# plt.xlabel("Comprimento [m]")
# plt.ylabel("Velocidade [m/s]")
# plt.title("Velocidade - Navier-Stokes")
# plt.grid()
# plt.figure()

# U = []
# for i in range(1, len(D)-1):
#     U.append(((P[i-1]-P[i])/dx) * ((D[i]**2) / (12 * mu)))

# plt.show()

x_p = np.linspace(0+dx, L-dx, num=nx, endpoint=True)
x_u = np.linspace(0, L, num=nu, endpoint=True)

mu = 1
dx = 1
w = get_D(D0)
Kfrac = np.zeros(len(w))
for i in range(0, len(w)):
    Kfrac[i] = (w[i]**2) / (12 * mu)

Kmeio_a = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_a[i] = (dx*2*1e-16) / mu

Kmeio_b = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_b[i] = (dx*2*1e-15) / mu
    
Kmeio_c = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_c[i] = (dx*2*1e-14) / mu
    
Kmeio_d = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_d[i] = (dx*2*1e-13) / mu
    
Kmeio_e = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_e[i] = (dx*2*1e-12) / mu
    
Kmeio_f = np.zeros(len(x_u))
for i in range(0, len(x_u)):
    Kmeio_f[i] = (dx*2*1e-11) / mu

# # plt.ylim(10e-30, 10e-7)
# # plt.xlim(-0.5, 10.5)
# plt.semilogy(x_u, Kfrac, label='Permeabilidade da Fratura')
# plt.semilogy(x_u, Kmeio_a, label='Permeabilidade do Meio (K) - 1e-16')
# plt.semilogy(x_u, Kmeio_b, label='K = 1e-15')
# plt.semilogy(x_u, Kmeio_c, label='K = 1e-14')
# plt.semilogy(x_u, Kmeio_d, label='K = 1e-13')
# plt.semilogy(x_u, Kmeio_e, label='K = 1e-12')
# plt.semilogy(x_u, Kmeio_f, label='K = 1e-11')
# plt.xlabel("Comprimento da Fratura [m]")
# plt.ylabel("Permeabilidade [m²]")
# plt.legend()
# plt.grid()
# # plt.savefig('semilogy.pdf', dpi=1200)
# plt.figure()

# plt.plot(x_u, Kfrac, label='Permeabilidade da Fratura')
# plt.plot(x_u, Kmeio_a, label='Permeabilidade do Meio (K) - 1e-16')
# plt.plot(x_u, Kmeio_b, label='K = 1e-15')
# plt.plot(x_u, Kmeio_c, label='K = 1e-14')
# plt.plot(x_u, Kmeio_d, label='K = 1e-13')
# plt.plot(x_u, Kmeio_e, label='K = 1e-12')
# plt.plot(x_u, Kmeio_f, label='K = 1e-11')
# plt.xlabel("Comprimento da Fratura [m]")
# plt.ylabel("Permeabilidade [m²]")
# plt.legend()
# plt.grid()
# # plt.savefig('plot.pdf', dpi=1200)
# plt.figure()

# plt.loglog(x_u, Kfrac, label='Permeabilidade da Fratura')
# plt.loglog(x_u, Kmeio_a, label='Permeabilidade do Meio (K) - 1e-16')
# plt.loglog(x_u, Kmeio_b, label='K = 1e-15')
# plt.loglog(x_u, Kmeio_c, label='K = 1e-14')
# plt.loglog(x_u, Kmeio_d, label='K = 1e-13')
# plt.loglog(x_u, Kmeio_e, label='K = 1e-12')
# plt.loglog(x_u, Kmeio_f, label='K = 1e-11')
# plt.xlabel("Comprimento da Fratura [m]")
# plt.ylabel("Permeabilidade [m²]")
# plt.legend()
# plt.grid()
# # plt.savefig('loglog.pdf', dpi=1200)
# plt.figure()

def get_normu():
    u_norm = np.zeros((len(results_u), len(results_u[0])))
    for j in range(0, len(results_u[0])):
        for i in range(0, len(results_p)):
            u_norm[i][j] = results_u[i][j] / results_u[0][j]
    return u_norm

u_norm = get_normu()
plt.plot(x_p, results_p)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Pressão [Pa]")
plt.grid()
plt.savefig('p_desac.pdf', dpi=1200)
plt.figure()

plt.plot(x_u, results_u)
# plt.plot(x_u, norm_u)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Velocidade [m/s]")
plt.grid()
plt.savefig('v_desac.pdf', dpi=1200)
plt.figure()

D = get_D(D0)
frac_perf = np.zeros(len(D))
for i in range(0, len(D)):
    frac_perf[i] = D[i]/2
    
plt.plot(x_u, frac_perf)
# plt.plot(x_u, norm_u)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Metade da Abertura [m]")
plt.grid()
plt.savefig('perfil_abertura.pdf', dpi=1200)
plt.figure()

plt.show()
np.savetxt("x_u.csv", x_u, delimiter=",")
np.savetxt("x_p.csv", x_p, delimiter=",")

with open('ns_u.csv', 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerows(results_u)

with open('ns_p.csv', 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerows(results_p)

print(datetime.datetime.now() - begin_time)