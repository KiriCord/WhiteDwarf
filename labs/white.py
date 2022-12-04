import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# среднее число электронов
Ye = 26 / 56
rho0 = 9.79e5 * (1 / Ye)  # г * (см ** -3)
m0 = 5.67e33 * (Ye ** 2)  # г
r0 = 7.72e8 * Ye  # см
# rho0 = 1


def gamma(x):
    return (x**2) / (3 * np.sqrt(1 + x ** 2))


def odes(r, y):
    rho, m = y
    #x = y[0] ** (1./3)

    drho_dr = - (m * rho) / (gamma(rho ** (1/3)) * r ** 2)
    dm_dr = (r ** 2) * rho

    return [drho_dr, dm_dr]


def whiteDwarf(rho_c, method='RK45'):
    solutions = []
    for i in range(len(rho_c)):
        rho_temp = rho_c[i]  # / rho0
        rSpan = [3e-14, 10]
        initCond = [rho_temp, 0]
        def rho_f(r, y): return y[0] - 5.13e-17
        rho_f.terminal = True
        sol = solve_ivp(odes, t_span=rSpan, y0=initCond,
                        method=method, events=rho_f)
        # sol.t - радиус
        # sol.y[0] - плотность
        # sol.y[1] - масса

        rhoValues = sol.y[0] * rho0
        mValues = sol.y[1] * m0
        rValues = sol.t * r0

        solutions.append((rValues[-1], mValues[-1], rhoValues[-1]))

    return solutions


def plotSolution(solutions):
    # for radius, mass, densities in solutions:
    #     plt.plot(mass, radius)

    plt.plot([solution[1] for solution in solutions], [solution[0]
             for solution in solutions])

    plt.axvline(solutions1[-1][1], linestyle="--",
                color="black", label="Предел Чандрасекара")
    plt.legend(fontsize=14)
    plt.ylabel("Радиус")
    plt.xlabel("Масса")
    plt.title("Зависимость масса—радиус для белых карликов")


# ran = [10e-3, 10e-2, 10e-1, 10e1, 10e2, 10e3, 10e4]
ran = np.logspace(-3, 5, num=30)

solutions1 = whiteDwarf(ran)

# print(len(solutions1[0][0]))
plotSolution(solutions1)
# print(solutions1)
# plt.plot(solutions1[1], solutions1[0])


plt.show()
