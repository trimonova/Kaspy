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

hx = 0.05
hy = 0.05

t_step = 0.005
T_exp = 10
Lx = 2
Ly = 2

Courant_number = t_step/k/hx**2 + t_step/k/hy**2
print(Courant_number)

N = int(Lx/hx) # количество ячеек вдоль оси х
M = int(Ly/hy)
#wells_with_Q = {(int(0.430/hx),int(0.430/hy)): -0.000003}
wells_with_Q = {}
#wells_with_Q = {(int(0.309/hx),int(0.309/hy)): -0.00001, (int(0.309/hx), int(0.551/hy)): -0.00001, (int(0.551/hx),int(0.309/hy)): -0.00001}
wells_coord = [(round(0.5/hx), round(1.1489/hy))]
P_well = [1450000]
print(wells_coord)

CP_dict = {}  # словарь, в котором ключами являются координаты точек с давлениями, а значения - значения этих давлений

for i in range(len(wells_coord)):
    CP_dict[wells_coord[i]] = P_well[i]

Pres = 1*10**5 # давление в пласте
Pbound = 1*10**5 #  давление на границе\

Pres_distrib = np.ones((N + 2, M + 2)) * Pres
P_add_hor = np.ones((N, 1)) * Pres
P_add_vert = np.ones((1, M + 2)) * Pres


def PorePressure_in_Time(Pres_distrib):

    A = np.zeros((N, N))
    B = np.zeros((N * M, 1))


    for n in range(1, N-1):
        A[n][n-1] = 1/hx**2
        A[n][n] = (-2/hx**2 - k/t_step - 2/hy**2)
        A[n][n+1] = 1/hx**2

    A[0][0] = (-2/hx**2 - k/t_step - 2/hy**2)
    A[0][1] = 1/hx**2
    A[N-1][N-1] = A[0][0]
    A[N-1][N-2] = A[0][1]

    A_sym = np.zeros((N, N))
    for n in range(0, N):
        A_sym[n][n] = 1/hy**2

    A_line_1 = np.hstack((A, A_sym, np.zeros((N, N * M - 3 * N)), A_sym))
    A_full = A_line_1.copy()

    for m in range(1, M - 1):
        A_line = np.hstack(
            (np.zeros((N, N * (m - 1))), A_sym, A, A_sym, np.zeros((N, N * M - (3 + (m - 1)) * N))))
        A_full = np.vstack((A_full, A_line))

    A_line_end = np.hstack((A_sym, np.zeros((N, N * M - 3 * N)), A_sym, A))
    A_full = np.vstack((A_full, A_line_end))

    j = 0
    for m in range(M):
        for n in range(N):
            if n == 0:
                # print(j, n, m)
                B[j][0] = -1/t_step * Pres_distrib[n][m]
            else:
                B[j][0] = -1/t_step * Pres_distrib[n][m]
            j += 1

    wells_coord.sort()
    wells_coord_reverse = wells_coord[:: -1]
    for well_coord_couple in wells_coord_reverse:
        A_well_column = A_full[:][(well_coord_couple[1] - 1) * N + well_coord_couple[0]]
        for cell_number in range(len(A_well_column)):
            if A_well_column[cell_number] != 0:
                B[cell_number] = B[cell_number] - A_well_column[cell_number] * CP_dict[well_coord_couple]
        A_full = np.delete(A_full, (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=0)
        A_full = np.delete(A_full,
                           (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=1)
        B = np.delete(B, (well_coord_couple[1] - 1) * N + well_coord_couple[0], axis=0)

    # print(np.shape(A_full), np.shape(B))
    P_new = np.linalg.solve(A_full, B)
    print(min(P_new), max(P_new))
    for well_coord_couple in wells_coord_reverse:
        P_new = np.insert(P_new, (well_coord_couple[1] - 1) * N + well_coord_couple[0],
                          CP_dict[well_coord_couple])
    # print(N_r, M_fi, N_r_oil, np.shape(P_new))
    P_new = P_new.reshape(N * M, 1)
    Pres_end = np.zeros((N, M))
    j = 0
    for m in range(M):
        for n in range(N):
            Pres_end[n][m] = P_new[j][0]
            j += 1

                    # for coord_key in wells_with_Q:
        #     if (n,m) == coord_key:
        #         B[n-1][0] = -V*beta/t_step*Pres_distrib[n][m]- alpha*coeff_1*(Pres_distrib[n][m-1] - 2*Pres_distrib[n][m] + Pres_distrib[n][m+1]) + wells_with_Q[coord_key]


    print(np.shape(Pres_end), np.shape(P_add_hor))
    #P_total = np.hstack((P_total,P_new))


    P_total = np.hstack((P_add_hor, Pres_end, P_add_hor))
    P_total = np.vstack((P_add_vert, P_total, P_add_vert))
    #print(np.shape(P_total))

    Pres_end = np.array(P_total.copy())

    return Pres_end


#----------------------------------------------------------------------------
for t in range(T_exp):
    print(t)
    Pres_distrib = PorePressure_in_Time(Pres_distrib)
    print(min(Pres_distrib.flat), max(Pres_distrib.flat))

X = np.zeros((N+2,M+2))
Y = np.zeros((N+2, M+2))
for m in range(M+2):
    for n in range(N+2):
        X[n][m] = n*hx
        Y[n][m] = m*hy

X_list = [i for i in X.flat]
Y_list = [j for j in Y.flat]
P_list = [k for k in Pres_distrib.flat]


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
