#  Можно задавать много точек с источниками и давлениями. Если задавать только давления (граничные условия) то задача устойчива при любых шагах времени и координаты,
# если задавать еще источники, то задача устойчива при каком-то соотношении t_step и hx, Q задается вроде бы в м3/с
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

perm = 2 * 10 ** (-15)  # м2 проницаемость
mu = 2 * 10 ** (-3)  # Па*с вязкость
fi = 0.2  # пористость
Cf = 10 ** (-9)  # сжимаемость флюида
Cr = 5 * 10 ** (-10)  # сжимаемость скелета
k = mu*fi*(Cf+Cr)/perm

alpha = 0.8*10**-12
beta = 0.17*10**-9
hx = 0.05
hy = 0.05

t_step = 0.005
T_exp = 5
Lx = 2
Ly = 2

Courant_number = t_step/k/hx**2 + t_step/k/hy**2
print(Courant_number)

N = int(Lx/hx) # количество ячеек вдоль оси х
M = int(Ly/hy)
#wells_with_Q = {(int(0.430/hx),int(0.430/hy)): -0.000003}
wells_with_Q = {}
#wells_with_Q = {(int(0.309/hx),int(0.309/hy)): -0.00001, (int(0.309/hx), int(0.551/hy)): -0.00001, (int(0.551/hx),int(0.309/hy)): -0.00001}
wells_with_P = {(round(0.5/hx), round(1.1489/hy)): 14.5*10**5}

Pres = 1*10**5 # давление в пласте
Pbound = 1*10**5 #  давление на границе\
P_bound_right = 5*10**5
P_bound_left = 1

def PorePressure_in_Time():
    Pres_distrib = np.ones((N + 2, M + 2)) * Pres  # пластовое давление во всей области на нулевом временном шаге
    indic = []

    for t in range(1, T_exp):
        P_add_hor = np.ones((N, 1)) * Pres
        P_add_vert = np.ones((1, M + 2)) * Pres
        P_total = np.ones((N, 1)) * Pres
        for m in range(1,M+1):
            A = np.zeros((N,N))
            B = np.zeros((N,1))

            for n in range(1, N-1):
                A[n][n-1] = 1/hx**2
                A[n][n] = (-2/hx**2 - k/t_step)
                A[n][n+1] = 1/hx**2

            A[0][0] = (-2/hx**2 - k/t_step)
            A[0][1] = 1/hx**2
            A[N-1][N-1] = A[0][0]
            A[N-1][N-2] = A[0][1]

            for n in range(1,N+1):
                if n == 1  or n == N:
                    B[n-1][0] = -k/t_step*Pres_distrib[n][m]- 1/hy**2*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1]) - 1/hx**2*Pbound

                else:
                    B[n-1][0] = -k/t_step*Pres_distrib[n][m]- 1/hy**2*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1])

                # for coord_key in wells_with_Q:
                #     if (n,m) == coord_key:
                #         B[n-1][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1]) + wells_with_Q[coord_key]

            for n in range(0, N):

                for coord_key in wells_with_P:
                    if (n,m-1) == coord_key:
                        indic.append(coord_key)
                    elif (n-1,m-1) == coord_key:
                        A[n][n - 1] = 0
                        B[n][0] = B[n][0] - 1/hx**2 * wells_with_P[coord_key]
                    elif (n+1,m-1) == coord_key:
                        A[n][n+1] = 0
                        B[n][0] = B[n][0] - 1/hx**2 * wells_with_P[coord_key]
            #print(type(indic))
            if indic != []:
                counter = 0
                for element in indic:
                    A = np.delete(A, element[0]-counter, axis=0)
                    A = np.delete(A, element[0]-counter, axis=1)
                    B = np.delete(B, element[0]-counter)
                    counter += 1


            P_new = np.linalg.solve(A,B)
    #
            if indic != []:
                counter = 0
                for element in indic:
                    P_new = np.insert(P_new,element[0]+counter,wells_with_P[element])
                    counter += 1
                indic = []
                P_new = P_new.reshape(N, 1)

            P_total = np.hstack((P_total,P_new))


        P_total = np.hstack((P_total, P_add_hor))
        P_total = np.vstack((P_add_vert, P_total, P_add_vert))

        Pres_distrib = np.array(P_total.copy())


    #---------------------------------------------------------------------------
        indic = []
        P_add_hor = np.ones((N + 2, 1)) * Pres
        P_add_vert = np.ones((1, M)) * Pres
        P_total = np.ones((1, M)) * Pres
        for n in range(1, N + 1):
            A = np.zeros((M, M))
            B = np.zeros((M, 1))
            for m in range(1, M - 1):
                A[m][m - 1] = 1/hy**2
                A[m][m] = (-2/hy**2 - k/t_step)
                A[m][m + 1] = 1/hy**2
            A[0][0] = (-2/hy**2 - k/t_step)
            A[0][1] = 1/hy**2
            A[M - 1][M - 1] = A[0][0]
            A[M - 1][M - 2] = A[0][1]

            for m in range(1, M + 1):
                if m == 1 or m == M:
                    B[m - 1][0] = -k/t_step*Pres_distrib[n][m]- 1/hx**2*(Pres_distrib[n-1][m] - 2*Pres_distrib[n][m] + Pres_distrib[n+1][m]) - 1/hy**2*Pbound

                else:
                    B[m - 1][0] = -k/t_step*Pres_distrib[n][m]- 1/hx**2*(Pres_distrib[n-1][m] - 2*Pres_distrib[n][m] + Pres_distrib[n+1][m])

                # for coord_key in wells_with_Q:
                #     if (n, m) == coord_key:
                #         B[m - 1][0] = -V * beta / t_step * Pres_distrib[n][m] - alpha * coeff_1 * (Pres_distrib[n][m - 1] - 2 * Pres_distrib[n][m] + Pres_distrib[n][m + 1]) + wells_with_Q[coord_key]


            for m in range(0, M):
                for coord_key in wells_with_P:
                    if (n-1,m) == coord_key:
                        indic.append(coord_key)

                    elif (n-1,m-1) == coord_key:
                        A[m][m - 1] = 0
                        B[m][0] = B[m][0] - 1/hy**2 * wells_with_P[coord_key]

                    elif (n-1,m+1) == coord_key:
                        A[m][m + 1] = 0
                        B[m][0] = B[m][0] - 1/hy**2 * wells_with_P[coord_key]

            if indic != []:
                counter = 0
                for element in indic:
                    A = np.delete(A, element[1]-counter, axis=0)
                    A = np.delete(A, element[1]-counter, axis=1)
                    B = np.delete(B, element[1]-counter)
                    counter += 1

            P_new = np.linalg.solve(A,B)

            if indic != []:
                counter = 0
                for element in indic:
                    P_new = np.insert(P_new, element[1]+counter, wells_with_P[element])
                    counter += 1
                indic = []
                P_new = P_new.reshape(M, 1)


            P_total = np.vstack((P_total, P_new.T))
        P_total = np.vstack((P_total, P_add_vert))
        P_total = np.hstack((P_add_hor, P_total, P_add_hor))
        Pres_distrib = np.array(P_total.copy())

    return P_total

#----------------------------------------------------------------------------

P_total = PorePressure_in_Time()
X = np.zeros((N+2,M+2))
Y = np.zeros((N+2, M+2))
for m in range(M+2):
    for n in range(N+2):
        X[n][m] = n*hx
        Y[n][m] = m*hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in P_total.flat]


CP_list = zip(X_list, Y_list, P_list)

if __name__ == '__main__':
    print(min(P_list), max(P_list))

    xi = np.linspace(min(X_list),max(X_list), 700)
    yi = np.linspace(min(Y_list), max(Y_list), 700)
    xig, yig = np.meshgrid(xi, yi)
    Pi = interpolate.griddata((X_list,Y_list), P_list, (xig, yig), method='cubic')

    levels = list(range(0,1500000,50000))
    fig = plt.figure()
    surf = plt.contourf(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi),linewidth=0.2, levels=levels)
    #ax = fig.gca(projection='3d')

    #surf = ax.plot_surface(xig, yig, Pi, cmap=cm.jet, antialiased=True, vmin=np.nanmin(Pi), vmax=np.nanmax(Pi), linewidth=0.2)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
