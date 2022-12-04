import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# среднее число электронов
Ye = 26 / 56
rho0 = 9.79 * 10e5 * (1 / Ye)  # г * (см ** -3)
m0 = 5.67 * 10e33 * (Ye ** 2)  # г
r0 = 7.72 * 10e8 * Ye  # см
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
    masses = []
    radi = []
    densities = []

    for i in range(len(rho_c)):
        rho_temp = rho_c[i]  # / rho0

        rSpan = [3e-10, 10]

        initCond = [rho_temp, 0]

        sol = solve_ivp(odes, t_span=rSpan, y0=initCond, method=method)

        # rhoValues = sol.y[0] * rho0
        # mValues = sol.y[1] * m0
        # rValues = sol.t * r0

        masses.append(sol.y[1][-1])
        densities.append(sol.y[0][-1])
        radi.append(sol.t[-1])

    return radi, masses, densities


def plotSolution(solutions):
    for radius, mass, densities in solutions:
        plt.plot(radius, mass)


ran = np.logspace(-3, 5, num=7)
solutions1 = whiteDwarf(ran)

# print(len(solutions1[0][0]))
# plotSolution(solutions1)
print(solutions1)
plt.plot(solutions1[1], solutions1[0])


plt.axvline(solutions1[1][-1], linestyle="--",
            color="black", label="Предел Чандрасекара")
plt.legend(fontsize=14)

plt.ylabel("Радиус")
plt.xlabel("Масса")
plt.title("Зависимость масса—радиус для белых карликов")
plt.show()
