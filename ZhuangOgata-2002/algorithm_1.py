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
    def __init__(self, archivo, M0=4.3, d=0.2, epsilon=10**(-3), max_iter=500):
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
        self.calculate_bandwidth()

        # 2. Set l = 0 and u^{(0)}(x,y) = 1
        self.u_xy = np.asarray([], dtype=float) 

        self.l = 0  # For iterations
        self.u_xy = np.append(self.u_xy, 1)

        condition = True
        while condition and self.l < self.max_iter:
            # 3. Using the maximum likelihood procedure, fit the conditional
            #   intensity function λ(t,x,y|Ht) = vu^{(1)}(x,y) + 
            #                               \sum_{k:tk<t}κ(Mk)g(t-tk)*f(x-xk,y-yk|Mk)
            # to the earthquake data.
            self.fit_conditional_intensity()

            # 4. Calculate ρj for each j=1,2,...,N
            temp_p = np.zeros([self.N])

            for j in range(self.N):
                lambda_j = self.evaluate_intensity_j(j, v, )
                temp_p[j] = self.calculate_pj(j, lambda_j)
            # 5. Calculate μ(x,y) and record as u^{l+1}(x,y)
            mu = self.calculate_mu_estim(self.X, self.Y, temp_p)
            self.u_xy = np.append(self.u_xy, mu)

            # 6. If max_{(x,y)}|u^{l+1}(x,y)-u^{l}(x,y)|> ε, where ε is a small positive
            # number , then set l = l + 1 and go to step 3. Otherwise, take 
            # v*u^{l+1}(x,y) as the background rate and stop.

            temp = {
                "iteration": self.l,
                "log L": self.log_likelihood,
                "v": self.v,
                "A": self.A,
                "c": self.c,
                "alpha": self.alpha,
                "p": self.p,
                "d": self.d}
            
            self.convergence_df.loc[len(self.convergence_df)] = temp

            self.difference = self.calculate_difference()
            if self.difference <= self.epsilon:
                condition = False
                break
            
            self.l +=1

        self.background_rate = self.u_xy[-1]
        print(self.convergence_df)
        
    def read_csv(self, file):
        return pd.read_csv(file)

    def prepare_data(self):
        self.df = self.df[self.df['Magnitude']>= self.M0].copy()

        # Tiempo
        self.df["Arrival_T"] = pd.to_datetime({
            "year": self.df["Year"],
            "month": self.df["Month"],
            "day": self.df["Day"],
            "hour": self.df["Hour"],
            "minute": self.df["Minute"],
            "second": self. df["Second"]
        })

        self.df = self.df.sort_values("Arrival_T").reset_index(drop=True)  # Ordenar

        self.T = (
            (self.df["Arrival_T"] - self.df["Arrival_T"].iloc[0])
            .dt.total_seconds()
        ).to_numpy()  # total time from any earthquake to the first earthquake

        self.df["Interarrival_T"] = self.df["Arrival_T"].diff().dt.total_seconds()
        self.T_inter = self.df["Interarrival_T"]
        self.T_total = self.T[-1]  # total time since the first earthquake

        self.N = len(self.df)
        self.M = self.df['Magnitude']
        self.X = self.df['Longitude']
        self.Y = self.df['Latitude']

    def calculate_bandwidth(self):
        # self.n_p is involved in the calculation of dj, but I set de degree value as the article
        self.dj = np.full(self.N, self.d)

    def kappa(self, M, A, alpha):
        return A * np.exp(alpha * (M - self.M0))

    def g(self, t, c, p):
        """
        t and c have to be in seconds
        """
        if t>0:
            return (p-1) * c**(p-1) * (t+c)**(-p)
        else:
            return 0

    def f(self, x, y, M, d, alpha):
        magnitude_factor = np.exp(alpha * (M - self.M0))

        denominator = 2 * np.pi * d * magnitude_factor

        numerator =  np.exp(- (x**2 + y**2) / (2 * d * magnitude_factor))

        return  numerator / denominator

    def log_likelihood(self, params):
        v, A, c, alpha, p = params
        log_history = 0

        for k in range(self.N):
            lambda_k = self.evaluate_intensity_j(k, v, A, c, alpha, p)
            log_history += np.log(lambda_k)

        # integral of the intensity
        integral_back = 1
        integral_offspring = 1

        integral = integral_back + integral_offspring

        l_L_lambda = log_history - integral
        return l_L_lambda

    def fit_conditional_intensity(self):
        """
        This function is for estimate the parameters
        """
        x0 = [self.v, self.A, self.c, self.alpha, self.p]

        result = minimize(lambda params: -self.log_likelihood(params), x0, method="Nelder-Mead")  # We must review other methods

        self.v = result.x[0]
        self.A = result.x[1]
        self.c = result.x[2]  # segundos
        self.alpha = result.x[3]
        self.p = result.x[4]

        self.log_likelihood_value = -result.fun

    def evaluate_intensity_j(self, j, v, A, c, alpha, p):
        background = v * self.u_xy[j]

        offspring = 0

        for i in range(j):

            Mi = self.M[i]

            delta_t = self.T[j] - self.T[i]
            delta_x = self.X[j] - self.X[i]
            delta_y = self.Y[j] - self.Y[i]

            offspring += (
                self.kappa(Mi, A, alpha)
                * self.g(delta_t, c, p)
                * self.f(
                    delta_x,
                    delta_y,
                    Mi,
                    self.dj[i],
                    alpha
                )
            )

        return background + offspring
    
    def calculate_pij(self, i, j, lambda_j):
        """
        proba of the jth eartquake being an offspring of ith event

        lambda_j : intensity function for the event j
                    it has to be calculated before
        """        
        Mi = self.M[i]
        di = self.d_j[i]
        delta_t = self.T[j] - self.T[i]
        delta_x = self.X[j] - self.X[i]
        delta_y = self.Y[j] - self.Y[i]

        pij = (self.kappa(Mi, self.A, self.alpha) 
               * self.g(delta_t, self.c, self.p) 
               * self.f(delta_x, delta_y, Mi, di, self.alpha)
               ) / lambda_j
        return pij
    
    def calculate_pj(self, j, lambda_j):
        """
        proba of the jth eartquake being an offspring in the process
        """
        pj = 0
        for i in range(j):  # from i=1 to j-1 in the article, but range goes from i=0 to j-1
            pij = self.calculate_pij(i,j, lambda_j)
            pj += pij

        return pj

    def calculate_mu_estim(self, x, y, p):
        temp = 0

        for j in range(self.N):
            delta_x = x - self.X[j]
            delta_y = y - self.Y[j]

            kdj = (2 * np.pi * self.d)**(-1) * np.exp(-(delta_x**2 + delta_y**2) / (2 * self.d**2))
            temp += (1-p) * kdj

        return temp / self.T_total

    def calculate_difference(self):
        return np.max(np.abs(self.u_xy[-1] - self.u_xy[-2]))


if __name__ =='__main__':
    earthquakes = "Earthquakes.csv"
    obj = ETAS_Declustering(earthquakes)
    print(obj.df.head())
    print(obj.convergence_df)
