import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve
import csv
import datetime
begin_time = datetime.datetime.now()

#Propriedades
rho = 1000 #kg/m³
cf =  0#5e-10 #1/Pa
mu = 1e-3 #Pa.s
K = 1e-15 #m²
Q = 0.000001 #m3/s
# Q = 1.25e-5

#Geometria
L = 10#Comprimento da Fratura
w0 = .001 #Abertura da entrada m

#Características de Tempo
dt = 1
tf = 1000000000000000
nt = int(tf / dt)

#Características de Malha
nx = 1000
nu = nx+1
dx = L/nx
ny = nx
dy = 0.1
tol = 1e-6

x_p = np.linspace(0+dx, L-dx, num=nx, endpoint=True)
x_u = np.linspace(0, L, num=nu, endpoint=True)

def get_x():
    x = np.zeros((nx+1, 1))
    for i in range(0, len(x)):
        x[i] = x[i-1] + dx
    return x

def get_w(w0):
    global w
    x = get_x()
    w = np.zeros(nx+1)
    for i in range(0, len(w)):
        w[i] = w0
    w = np.linspace(w0, w0/10, num = int(len(w)), endpoint=True)
    # for i in range(0, int(len(x)/1.666)):
    #     w[i] = (w0*(1-((x[i]/(L/1.666))**2)))
    # for i in range(int(len(x)/1.666), len(x)):
    #     w[i] = 0
    # for i in range(0, int(len(x)/1)):
    #     w[i] = (w0*(1-((x[i]/(L/1))**2)))
    return w
    
def get_wp():
    global wp
    w = get_w(w0)
    wp = np.zeros(nx)
    for i in range(0, len(wp)):
        wp[i] = ((w[i]+w[i+1])/2)
    return wp

def get_dy():
    wp = get_wp()
    dy = np.zeros(nx)
    dy[0] = 0.1
    for i in range(0, len(dy)-1):
        dy[i+1] = dy[0]
        # dy[i+1] = dy[0] + (wp[0] - wp[i+1])/2
    return dy

#Cálculo de Coeficientes
def get_Ae():
    w = get_w(w0)
    Ae = np.zeros(nx)
    for i in range(0, nx):
        Ae[i] =  (((w[i+1]**3)) / (12 * mu * dx)) #* rho
    return Ae

def get_Aw():
    w = get_w(w0)
    Aw = np.zeros(nx)
    for i in range(0, nx):
        Aw[i] = (((w[i]**3)) / (12 * mu * dx)) #* rho
    return Aw

def get_Ap0():
    wp = get_wp()
    Ap0 = np.zeros(nx)
    for i in range(0, nx):
        Ap0[i] = ((dx * wp[i] * cf) / dt)
    return Ap0

def get_An():
    An = np.zeros(nx)
    dy = get_dy()
    for i in range(0, ny):
        An[i] = ((K * dx) / (mu * dy[i])) #* rho
    return An

def get_As():
    As = np.zeros(nx)
    dy = get_dy()
    for i in range(0, ny):
        As[i] = ((K * dx) / (mu * dy[i])) #* rho
    return As

def get_Ap():
    Ap = np.zeros(nx)
    Ae = get_Ae()
    Aw = get_Aw()
    An = get_An()
    As = get_As()
    Ap0 = get_Ap0()
    for i in range(0, nx):
        Ap[i] = Ae[i] + Aw[i] + An[i] + As[i] + Ap0[i]
    return Ap

#Construção da Matriz de Coeficientes
def get_A():
    A = np.zeros((nx, nx))
    Ap = get_Ap()
    Ae = get_Ae()
    Aw = get_Aw()
    An = get_An()
    As = get_As()
    Ap0 = get_Ap0()
    for i in range(0, nx):
         A[i,i] = Ap[i]
    for j in range(1, nx):
         A[j, j-1] = -Aw[j]
    for k in range(0, nx-1):
         A[k, k+1] = -Ae[k]
    # A[0, 0] = Ae[0] + An[0] + As[0] + Ap0[0]
    A[-1, -1] = Aw[-1] + An[-1] + As[-1] + Ap0[-1]
    A[0, 0] = Ap0[0] + Ae[0] + 2*Aw[0] + An[0] + As[0]
    # A[-1, -1] = Ap0[-1] + Aw[-1] + 2*Ae[-1] + An[-1] + As[-1]
    return A

#Construção do Termo Fonte
def get_Pn():
    Pn = np.zeros((ny,1))
    for i in range(0, len(Pn)):
        Pn[i] = 0
    return Pn

def get_Ps():
    Ps = np.zeros((ny,1))
    for i in range(0, len(Ps)):
        Ps[i] = 0
    return Ps

def get_b(Pf0, b0):
    global Pbw
    b = b0.copy()
    Ae = get_Ae()
    Aw = get_Aw()
    Ap0 = get_Ap0()
    Pn = get_Pn()
    Ps = get_Ps()
    An = get_An()
    As = get_As()
    Pbw = 5e6
    Pbe = 0
    for i in range(0, nx):
        b[i] = Pn[i]*An[i] + Ps[i]*As[i] + Ap0[i]*Pf0[i]
    # b[0] =  Pn[0]*An[0] + Ps[0]*As[0] + Q + Ap0[0]*Pf0[0]
    b[-1] =  Pn[-1]*An[-1] + Ps[-1]*As[-1] + 0 + Ap0[-1]*Pf0[-1]
    b[0] = Ap0[0]*Pf0[0] + Pn[0]*An[0] + Ps[0]*As[0] + 2*Aw[0]*Pbw
    # b[-1] = Ap0[-1]*Pf0[-1] + Pn[-1]*An[-1] + Ps[-1]*As[-1] + 2*Ae[-1]*Pbe
    return b

#Solver
def solve():
    Pf0 = np.full((nx,1), 0)
    for i in range(0, nx):
        Pf0[i] = 0
    b0 = np.zeros((nx,1))
    A = get_A()
    t = 0
    while (t <= tf):
        b = get_b(Pf0, b0)
        lu, piv = lu_factor(A)
        Pf = lu_solve((lu, piv), b)
        # Pf = np.linalg.solve(A, b)
        if (np.linalg.norm(Pf-Pf0,np.inf) < tol):
            print(t)
            break
        else:
            Pf0 = Pf.copy()
            t = t + dt
    if np.allclose(np.dot(A, Pf), b) == True:
        print('Solução Encontrada')
    else:
        print('Não foi obtida uma solução')
    return Pf

def get_Plist():
    global K, Klist
    Plist = []
    Klist = []
    K = 1e-15
    while K <= 1e-15:
        Klist.append(K)
        Pf = solve()
        Plist.append(Pf)
        K = K*10
    return Plist

#Análise de inclinação das curvas para valores variados de permeabilidade
def get_m(Pf):
    x = get_x()
    a0 = round((1/2)*L)
    n = np.abs(x - a0).argmin()
    m = []
    for i in range(0, n):
        m.append((Pf[i+1] - Pf[i])/(x[i+1] - x[i]))
    return m

def checkCurve(m):
    check = []
    for i in range(0, len(m)-1):
        check.append(m[i] < m[i+1])
    return check[0]

def get_grad():
    Ae = get_Ae()
    w = get_w(w0)
    grad = []
    for i in range(0, len(Ae)-1):
        grad.append((w[i] - w[i+1])/((dx)))
    gradm = np.mean(grad)
    return gradm

def get_Curve():
    global w0
    global K
    gradm = []
    K = 1e-16
    Kv = []
    wv = []
    while K <= 1e-11:
        w0 = 0.2
        while w0 > 0:
            Pf = solve()
            m = get_m(Pf)
            if checkCurve(m) == True:
                gradm.append(get_grad())
                Kv.append(K)
                wv.append(w0)
                # print(K)
                # print(gradm[-1])
                break
            else:
               w0 = (99*w0)/100
            # print(w0)
        # print(K)
        K = K*1.5
    return gradm, Kv, wv

def write_txt():
    gradm, Kv, wv = get_Curve()
    with open("grad.txt", "ab") as f:
        np.savetxt(f, gradm)
    with open("K.txt", "ab") as f:
        np.savetxt(f, Kv)
    with open("w.txt", "ab") as f:
        np.savetxt(f, wv)
    return

# gradm, Kv, wv = get_Curve()

# gradm = []
# Kv = []
# wv = []
# for line in open('grad.txt', 'r'):
#   values = [float(s) for s in line.split()]
#   gradm.append(values[0])

# for line in open('K.txt', 'r'):
#   values = [float(s) for s in line.split()]
#   Kv.append(values[0])
  
# for line in open('w.txt', 'r'):
#   values = [float(s) for s in line.split()]
#   wv.append(values[0])
  
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.plot(gradm, Kv)
# plt.xlabel("dw/dx")
# plt.ylabel("K")
# plt.figure()
# plt.semilogy(gradm, Kv)
# plt.xlabel("dw/dx")
# plt.ylabel("K")
# plt.figure()
# plt.semilogx(gradm, Kv)
# plt.xlabel("dw/dx")
# plt.ylabel("K")
# plt.figure()
# plt.loglog(gradm, Kv)
# plt.xlabel("dw/dx")
# plt.ylabel("K")
# plt.figure()

# plt.semilogy(wv, Kv)
# plt.xlabel("w")
# plt.ylabel("K")
# plt.figure()

# w = get_w(w0)
# x = np.zeros((nx,1))
# for i in range(0, len(x)):
#     x[i] = x[i-1] + dx
# Pf = solve()
    
# plt.ticklabel_format(axis="y", useOffset=False, style="sci")
# plt.plot(x, Pf)
# plt.xlabel("Comprimento da Fratura [m]")
# plt.ylabel("Pressão [Pa]")
# plt.title("Pressão - Poiseuille")
# plt.grid()
# plt.figure()

# def get_u(Pf):
#     w = get_w(w0)
#     u = np.zeros((len(w),1))
#     for i in range(1, len(w)-1):
#         u[i] = (((Pf[i-1]-Pf[i])/dx) * ((w[i]**2) / (12 * mu)))
#     u[0] = Q/w0
#     return u

def get_u(Pf):
    w = get_w(w0)
    u = np.zeros((len(w),1))
    for i in range(1, len(w)-1):
        u[i] = (((Pf[i-1]-Pf[i])/dx) * ((w[i]**2) / (12 * mu)))
    u[0] = (((Pbw-Pf[0])/dx) * ((w[0]**2) / (12 * mu)))
    return u

# def get_ulist(Plist):
#     ulist = []
#     u = np.zeros((len(w),1))
#     i = 0
#     while i < int(len(Plist)):
#         u = np.zeros((len(w),1))
#         for j in range(1, len(w)-1):
#             u[j] = (((Plist[i][j-1]-Plist[i][j])/dx) * ((w[j]**2) / (12 * mu)))
#             u[0] = Q/w0
#         ulist.append(u)
#         i = i + 1
#     return ulist

def get_ulist(Plist):
    w = get_w(w0)
    ulist = []
    u = np.zeros((len(w),1))
    i = 0
    while i < int(len(Plist)):
        u = np.zeros((len(w),1))
        for j in range(1, len(w)-1):
            u[j] = (((Plist[i][j-1]-Plist[i][j])/dx) * ((w[j]**2) / (12 * mu)))
            u[0] = (((Pbw-Plist[i][0])/(dx/2)) * ((w[0]**2) / (12 * mu)))
        ulist.append(u)
        i = i + 1
    return ulist

# Pf = solve()
# u = get_u(Pf)
# plt.plot(u)
# plt.xlabel("Comprimento da Fratura [m]")
# plt.ylabel("Velocidade [m/s]")
# plt.title("Velocidade - Poiseuille")
# plt.grid()
# plt.figure()

Plist = get_Plist()
ulist = get_ulist(Plist)

# U = -((w0**2)/(12*mu))*((-Pbw)/(L))

# def get_m_ponto(u):
#     w = get_w(w0)
#     m_ponto = np.zeros(len(u))
#     for i in range(0, len(u)):
#         m_ponto[i] = (w[i]) * u[i]
#         m_ponto[i] = round(m_ponto[i], 6)
#     return m_ponto

# m_ponto = get_m_ponto(u)

# plt.plot(m_ponto)
# plt.grid()
# plt.figure()

# x = get_x()
# plt.plot(x, w)
# plt.figure()

def get_leakoff():
    Pf = solve()
    Pn = get_Pn()
    Ps = get_Ps()
    leakoff_n = np.zeros((nx,1))
    leakoff_s = np.zeros((nx,1))
    for i in range(0, nx):
        leakoff_n[i] = ((K / mu) * ((Pf[i] - Pn[i]) / dy)) * dx
        leakoff_s[i] = ((K / mu) * ((Pf[i] - Ps[i]) / dy)) * dx
    leakoff = np.zeros((nx,1))
    for j in range(0, nx):
        leakoff[j] = leakoff_n[j] + leakoff_s[j]
    return leakoff

leakoff = get_leakoff()
# leakoff_total = np.sum(leakoff)
# diff_leakoff_velocidade = (Q - leakoff_total)
# print(diff_leakoff_velocidade)
# plt.plot(x, leakoff)
# plt.grid()
# plt.figure()

# plt.show()

# for i in range(0, len(m_ponto)):
#     print(m_ponto[i] - leakoff[i])

Klist = np.array(Klist)

def get_results_u(ulist):
    results_u = np.zeros((len(ulist[0]), len(ulist)))
    for h in range(0, len(ulist)):
        for i in range(0, len(ulist[0])):
            results_u[i][h] = ulist[h][i][0]
    return results_u

def get_results_p(Plist):
    results_p = np.zeros((len(Plist[0]), len(Plist)))
    for h in range(0, len(Plist)):
        for i in range(0, len(Plist[0])):
            results_p[i][h] = Plist[h][i][0]
    return results_p

results_u = get_results_u(ulist)
results_p = get_results_p(Plist)

def get_flow(results_u, w):
    flow = np.zeros((len(results_u), len(results_u[0])))
    i = 0
    while i < int(len(results_u[0])):
        for j in range(0, len(results_u)):
            flow[j][i] = results_u[j][i] * w[j] * 1000
        i = i + 1
    return flow

# flow = get_flow(results_u, w)
# plt.plot(flow)
# plt.figure()


plt.plot(x_p, results_p)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Pressão [Pa]")
plt.title("Abertura Variável")
plt.grid()
# plt.savefig('p_nx3_var.pdf', dpi=1200)
plt.figure()

plt.ticklabel_format(axis="y", useOffset=True, style="sci")
plt.plot(x_u, results_u)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Velocidade [m/s]")
plt.title("Abertura Variável nx=1000")
plt.grid()
# plt.savefig('v_nx1000_var.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.show()

np.savetxt("K.csv", Klist, delimiter=",")

with open('poiseuille_u.csv', 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerows(results_u)
        
with open('poiseuille_p.csv', 'w') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerows(results_p)

print(datetime.datetime.now() - begin_time)