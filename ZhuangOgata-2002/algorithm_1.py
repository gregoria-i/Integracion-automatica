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
    def __init__(self, archivo, M0=4.3, d=0.02, epsilon=1e-3, max_iter=20):  # change M0=4.3
        self.archivo = archivo
        self.M0 = M0
        self.d = d
        self.epsilon = epsilon  
        self.N = 0
        self.max_iter = max_iter
        self.convergence_df = pd.DataFrame(columns = ["iteration", "log L", "v", "A", "c", "alpha", "p", "d"])

        np.random.seed(121)

        self.v = np.random.uniform(0.1, 1.0)  
        self.A = np.random.uniform(0.1, 1.0)  
        self.c = np.random.uniform(0.1, 1.0)  
        self.alpha = np.random.uniform(0.1, 1.0)  
        self.p = np.random.uniform(1.1, 2.0)  

        self.df = self.read_csv(self.archivo)

        self.prepare_data()
        
        # 1. Given a preliminary parameter np, say 20, calculate the bandwidth dj
        #   for each event (tj, xj, yj, Mj) 
        self.n_p = 20  # at least np other earthquakes

        self.calculate_bandwidth()

        # 2. Set l = 0 and u^{(0)}(x,y) = 1
        self.u_xy = np.ones(self.N)  # for each one of the earthquakes
        self.u_xy_new = np.ones(self.N)  # with the same size 

        self.l = 0  # For iterations

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
                #print(f"j:{j}, v:{self.v}, A:{self.A}, c:{self.c}, alpha:{self.alpha}, p:{self.p}")
                lambda_j = self.evaluate_intensity_j(j, self.v, self.A, self.c, self.alpha, self.p)  # lambdaj has to be >0
                temp_p[j] = self.calculate_pj(j, lambda_j)  # We have N lambdaj

            # 5. Calculate μ(x,y) and record as u^{l+1}(x,y)
            mu = self.calculate_mu_estim(self.X, self.Y, temp_p)  # len(mu) = self.N
            self.u_xy_new = mu

            # 6. If max_{(x,y)}|u^{l+1}(x,y)-u^{l}(x,y)|> ε, where ε is a small positive
            # number , then set l = l + 1 and go to step 3. Otherwise, take 
            # v*u^{l+1}(x,y) as the background rate and stop.
            
            params = [self.v, self.A, self.c, self.alpha, self.p]
            log_L = self.log_likelihood(params)

            temp = {
                "iteration": self.l,
                "log L": log_L,
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
            self.u_xy = self.u_xy_new.copy()

        self.save_results()
        
    def read_csv(self, file):
        return pd.read_csv(file)

    def prepare_data(self):
        self.df = self.df[self.df['Magnitude']>= self.M0].copy()

        # Time
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
            .dt.total_seconds()/ 86400 # days
        ).to_numpy()  # total time from any earthquake to the first earthquake

        self.df["Interarrival_T"] = self.df["Arrival_T"].diff().dt.total_seconds() / 86400
        self.T_inter = self.df["Interarrival_T"]
        self.T_total = self.T[-1]  # total time since the first earthquake

        self.N = len(self.df)
        self.M = self.df['Magnitude']
        self.X = self.df['Longitude']
        self.Y = self.df['Latitude']

    def calculate_bandwidth(self):
        """
        self.n_p is involved in the calculation of dj, but I set de degree value as the article
        """
        self.dj = np.full(self.N, self.d)

    def fit_conditional_intensity(self):
        """
        This function is for estimate the parameters and updates the global parameters
        """
        x0 = [self.v, self.A, self.c, self.alpha, self.p]

        epsilon = 1e-3
        bounds = [
            (0 + epsilon, None),  # v
            (0 + epsilon, 1 - epsilon),  # A
            (0 + epsilon, None),  # c
            (None, None),  # alpha
            (1 + epsilon, None),  # p
        ]

        def neg(x0):
            return -self.log_likelihood(x0)

        result = minimize(neg, x0, method="Nelder-Mead", bounds=bounds, tol=1e-3)

        self.v = result.x[0]
        self.A = result.x[1]
        self.c = result.x[2]  # days
        self.alpha = result.x[3]
        self.p = result.x[4]

    def log_likelihood(self, params):
        """
        the internal functions are evaluated with params, not with the global variables
        """
        v, A, c, alpha, p = params
        log_history = 0

        for k in range(self.N):
            lambda_k = self.evaluate_intensity_j(k, v, A, c, alpha, p)
            log_history += np.log(lambda_k)

        integral_back = v * self.T_total  # integral of the intensity

        integral_offspring = 0.0  # Offspring

        for i in range(self.N):

            Mi = self.M[i]

            T_remaining = self.T_total - self.T[i]  # Time available after event i

            integral_g = (1 - (c / (T_remaining + c)) ** (p - 1)) # Integral of g(t - ti) from ti to T

            integral_offspring += (self.kappa(Mi, A, alpha)* integral_g)


        integral = integral_back + integral_offspring

        l_L_lambda = log_history - integral  # l_L_lambda its a numpy.float64
        return l_L_lambda

    def evaluate_intensity_j(self, j, v, A, c, alpha, p):
        """
        this function has to return lambda()> 0
        """

        background = v * self.u_xy[j]  # >0

        offspring = 0

        for i in range(j):

            Mi = self.M[i]
            delta_t = self.T[j] - self.T[i]
            delta_x = self.X[j] - self.X[i]
            delta_y = self.Y[j] - self.Y[i]

            k_evaluated = self.kappa(Mi, A, alpha)
            g_evaluated = self.g(delta_t, c, p)  # >0
            f_evaluated = self.f(delta_x, delta_y, Mi, self.dj[i], alpha)  # sometimes =0

            offspring += k_evaluated * g_evaluated * f_evaluated  

        result = background + offspring
        return result

    def kappa(self, M, A, alpha):
        return A * np.exp(alpha * (M - self.M0))

    def g(self, t, c, p):
        """
        t and c have to be in days
        """
        if (t>0 and p > 1 and c>0):
            return (p-1) * c**(p-1) / (t+c)**(p)
        else:
            return 0

    def f(self, x, y, M, d, alpha):
        magnitude_factor = np.exp(alpha * (M - self.M0))  # >0

        denominator = 2 * np.pi * d * magnitude_factor  # >0

        numerator =  np.exp(- (x**2 + y**2) / (2 * d * magnitude_factor))  # >0

        return  numerator / denominator
    
    def calculate_pij(self, i, j, lambda_j):
        """
        proba of the jth eartquake being an offspring of ith event

        lambda_j : intensity function for the event j
                    it has to be calculated before
        """        
        Mi = self.M[i]
        di = self.dj[i]
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
        return np.max(np.abs(self.u_xy_new - self.u_xy))

    def evaluate_u_over_grid(self, X, Y):
        p = self.p
        Z = self.calculate_mu_estim(X, Y, p)
        return Z
    
    def save_results(self):
        self.convergence_df.to_csv("Convergence_table.csv", index=False)
        self.u_xy.name = "u_xy"
        self.u_xy.to_csv("U_xy.csv", index=False)
        

if __name__ =='__main__':
    earthquakes = "Earthquakes.csv"
    obj = ETAS_Declustering(earthquakes)
