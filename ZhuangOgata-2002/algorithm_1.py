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


class ETAS_Declustering:
    def __init__(self, archivo, M0=4.3, d=0.2, epsilon=10**(-3), max_iter=100):
        self.archivo = archivo
        self.M0 = M0
        self.d = d
        self.epsilon = epsilon  
        self.N = 0
        self.max_iter = max_iter
        self.convergence_df = pd.DataFrame(columns = ["iteration", "log L", "v", "A", "c", "alpha", "p", "d"])
        self.background_rate = 1

        self.df = self.read_csv(self.archivo)

        self.prepare_data()
        
        # 1. Given a preliminary parameter np, say 20, calculate the bandwidth dj
        #   for each event (tj, xj, yj, Mj) 
        self.n_p = 20  # at least np other earthquakes
        d = self.calculate_bandwidth()

        # 2. Set l = 0 and u^{(0)}(x,y) = 1
        self.u_xy = np.asarray([], dtype=int)  # Here we have to check if u is an array of int

        self.l = 0  # For iterations
        self.u_xy = np.append(self.u_xy, 1)

        condition = True
        while condition == True and self.l < self.max_iter:
            # 3. Using the maximum likelihood procedure, fit the conditional
            #   intensity function λ(t,x,y|Ht) = vu^{(1)}(x,y) + 
            #                               \sum_{k:tk<t}κ(Mk)g(t-tk)*f(x-xk,y-yk|Mk)
            # to the earthquake data.
            self.fit_conditional_intensity()

            # 4. Calculate ρj for each j=1,2,...,N
            self.calculate_pj()

            # 5. Calculate μ(x,y) and record as u^{l+1}(x,y)
            mu = self.calculate_mu_estim()
            self.u_xy.append(self.u_xy, mu)

            # 6. If max_{(x,y)}|u^{l+1}(x,y)-u^{l}(x,y)|> ε, where ε is a small positive
            # number , then set l = l + 1 and go to step 3. Otherwise, take 
            # v*u^{l+1}(x,y) as the background rate and stop.

            self.difference = self.calculate_difference()

            self.temp = {
                "iteration": self.l,
                "log L": self.log_likelihood,
                "v": self.v,
                "A": self.A,
                "c": self.c,
                "alpha": self.alpha,
                "p": self.p,
                "d": self.d}
            self.convergence_df[self.l] = self.temp

            if self.difference <= self.epsilon:
                condition = False
                break
            
            self.l +=1

        self.background_rate = self.u_xy[-1]
        print(self.convergence_df)

    def read_csv(self, file):
        return pd.read_csv(file)

    def prepare_data(self):
        self.df = self.df[self.df['Magnitude']>= self.M0]
        self.N = len(self.df)
        self.M = self.df['Magnitude']
        self.X = self.df['Longitude']
        self.Y = self.df['Latitude']

        # Tiempo
        self.df["Arrival_T"] = pd.to_datetime({
            "year": self.df["Year"],
            "month": self.df["Month"],
            "day": self.df["Day"],
            "hour": self.df["Hour"],
            "minute": self.df["Minute"],
            "second": self. df["Second"]
        })
        self.T = self.df["Arrival_T"]

        self.df["Interarrival_T"] = self.df["Arrival_T"].diff()
        self.T_total = self.df['Arrival_T'].iloc[-1] - self.df['Arrival_T'].iloc[0]

    def calculate_bandwidth(self):
        # self.n_p is involved in the calculation of dj, but I set de degree value as the article
        self.dj = np.full(self.N, self.d)

    def kappa(self):
        return A * np.exp(alpha * (M - self.M0))

    def g(self, t, c, p):
        if t>0:
            return (p-1) * c**(p-1) * (t+c)**(-p)
        else:
            return 0

    def f(self):
        gau = (2 * np.pi * d * np.exp(alpha * (M - self.M0)))**(-1) * np.exp(- (x**2 + y**2) / (2 * d * np.exp(alpha * (M - self.M0))))
        return  gau

    def fit_conditional_intensity(self, t, x, y, Ht, u, idx, v, ):
        main_shocks_intensity = v * u[idx]  # this is the mu estimator
        other_shocks_intensity = 0
        for tk in Ht:
            if tk < t:
                other_shocks_intensity += self.kapa(M[k]) * g(t-tk) * f(x-x[k], y-y[k], M[k], d, self.M0, alpha)

        lamb = main_shocks_intensity + other_shocks_intensity

        l_lamb = np.log(lamb)
        x0 = [0]
        neg_func = -l_lamb
        result = minimize(neg_func, x0, method='BFGS')
        return result.x[0]

    def calculate_pj(self):
        pij = np.array([N])
        pij = self.kappa(Mi) * self.g(tj-ti) * self.f(xj-xi, yj-yi, Mi) / fci(tj, xj, yj, H_tj)
        pj = np.sum(pij) 
        return pj

    def calculate_mu_estim(self, x, y):
        kdj = (2 * np.pi * d)**(-1) * np.exp(-(x**2 + y**2) * (2 * d**2)**(-1))
        numerador = np.sum((1-pj) * kdj)
        mu_estim = T**(-1)  * numerador

        return mu_estim

    def calculate_difference(self):
        pass



if __name__ =='__main__':
    earthquakes = "Earthquakes.csv"
    obj = ETAS_Declustering(earthquakes)
    print(obj.df.head())
    print(obj.convergence_df)