"""
algorithm_2.py
@date: Aug-15-2026

With the estimated background intensity and the branching structure from
algorithm 1, Zhuang-Ogata proposed the following declustering algorithm
"""
import numpy as np

# 1. For each earthquake j=1, 2, ..., N, calculate probability ρj from
#   the final solution in algorithm 1.
N = 10
p = []
for j in range(0, N-1):
    p[j] = 0
# 2. Generate N uniform random numbers U1, U2, ...UN in [0,1].
values_Ux = np.random.uniform(0, 1, size=N)

# 3. If Uj < 1- ρj, then keep the jth event; otherwise, delete it from
#   the catalog as an offspring.
for j in range(0, N-1):
    if values_Ux[j] < 1 -p[j]:
        pass
    else: 
        values_Ux.pop(j)
