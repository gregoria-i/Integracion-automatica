"""
algorithm_1.py
@date: Aug-15-2026

To estimate the background rate from seismic data, Zhuang & Ogata use the
following iterative algorithm that simultaneously estimates the
background rate and the branching structure.
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize


def read_csv(file):
    return pd.read_csv(file)

def calculate_bandwidth(df):  # , center, n_p):
    d = [1] * len(df)   # Asignado provisionalmente por decisión con mi asesor
    return d

def kapa(M, A, alpha, M0):
    return A * np.exp(alpha * (M - M0))

def g(t, c, p):
    if t>0:
        return (p-1) * c**(p-1) * (t+c)**(-p)
    else:
        return 0

def f(x, y, M, d, M0, alpha):
    gau = (2 * np.pi * d * np.exp(alpha * (M - M0)))**(-1) * np.exp(- (x**2 + y**2) / (2 * d * np.exp(alpha * (M - M0))))
    return  gau

def fit_conditional_intensity(t, x, y, Ht, u, idx, v, ):
    main_shocks_intensity = v * u[idx]  # this is the mu estimator
    other_shocks_intensity = 0
    for tk in Ht:
        if tk < t:
            other_shocks_intensity += kapa(M[k]) * g(t-tk) * f(x-x[k], y-y[k], M[k], d, M0, alpha)

    lamb = main_shocks_intensity + other_shocks_intensity

    l_lamb = np.log(lamb)
    x0 = [0]
    neg_func = -l_lamb
    result = minimize(neg_func, x0, method='BFGS')
    return result.x[0]

def calculate_pj():
    pij = np.array([N])
    pij = kapa(Mi) * g(tj-ti) * f(xj-xi, yj-yi, Mi) / fci(tj, xj, yj, H_tj)
    pj = np.sum(pij) 
    return pj

def calculate_mu_estim(x, y):
    kdj = (2 * np.pi * d)**(-1) * np.exp(-(x**2, y**2) * (2d**2)**(-1))
    numerador = np.sum((1-pj) * kdj)
    mu_estim = T**(-1)  * numerador

    return mu_estim


if __name__ =='__main__':
    earthquakes = "Earthquakes.csv"
    df = read_csv(earthquakes)

    # 1. Given a preliminary parameter np, say 20, calculate the bandwidth dj
    #   for each event (tj, xj, yj, Mj) 
    n_p = 20  # at least np other earthquakes

    d = calculate_bandwidth(df)  #, center=(0,0), n_p=n_p)

    # 2. Set l = 0 and u^{(0)}(x,y) = 1
    u_xy = []

    l = 0
    u_xy[0] = 1

    N=20
    p = []

    # 3. Using the maximum likelihood procedure, fit the conditional
    #   intensity function λ(t,x,y|Ht) = vu^{(1)}(x,y) + 
    #                               \sum_{k:tk<t}κ(Mk)g(t-tk)*f(x-xk,y-yk|Mk)
    # to the earthquake data.
    while condition == True:
        fci = fit_conditional_intensity(t, x, y, Ht, v=1)

        # 4. Calculate ρj for each j=1,2,...,N

        for j in range(1,N):
            p[j] = calculate_pj()


        # 5. Calculate μ(x,y) and record as u^{l+1}(x,y)
            u_xy[l+1] = calculate_mu_estim()

        # 6. If max_{(x,y)}|u^{l+1}(x,y)-u^{l}(x,y)|> ε, where ε is a small positive
        # number , then set l = l + 1 and go to step 3. Otherwise, take 
        # v*u^{l+1}(x,y) as the background rate and stop.
            epsilon = 10**(-3)
            v= 20
            if np.abs(u_xy[l+1]-u_xy[l])>epsilon:
                l = l+1

            else:
                v * u_xy[l+1]
                condition = False

