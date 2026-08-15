"""
algorithm_1.py
@date: Aug-15-2026

To estimate the background rate from seismic data, Zhuang & Ogata use the
following iterative algorithm that simultaneously estimates the
background rate and the branching structure.
"""
import numpy as np

# 1. Given a preliminary parameter np, say 20, calculate the bandwidth dj
#   for each event (tj, xj, yj, Mj) 
N=20
n_p = 20
p = []
def calculate_bandwidth():
    pass 

for event_j in range(1, N):
    dj = calculate_bandwidth(event_j)

# 2. Set l = 0 and u^{(0)}(x,y) = 1
l = 0
u_xy = []

u_xy[0] = 1

# 3. Using the maximum likelihood procedure, fit the conditional
#   intensity function λ(t,x,y|Ht) = vu^{(1)}(x,y) + 
#                               \sum_{k:tk<t}κ(Mk)g(t-tk)*f(x-xk,y-yk|Mk)
# to the earthquake data.

# 4. Calculate ρj for each j=1,2,...,N
def calculate_pj():
    pass

for j in range(1,N):
    p[j] = calculate_pj()

# 5. Calculate μ(x,y) and record as u^{l+1}(x,y)

def calculate_mu():
    pass

u_xy[l+1] = calculate_mu()

# 6. If max_{(x,y)}|u^{l+1}(x,y)-u^{l}(x,y)|> ε, where ε is a small positive
# number , then set l = l + 1 and go to step 3. Otherwise, take 
# v*u^{l+1}(x,y) as the background rate and stop.
epsilon = 10^-3
v= 20
if np.abs(u_xy[l+1]-u_xy[l])>epsilon:
    l = l+1
else:
    v * u_xy[l+1]
