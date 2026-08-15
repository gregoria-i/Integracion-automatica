"""
algorithm_3.py
@date: Aug-15-2026

Zhuang & Ogata use the following algorithm to provide a stochastic
realization of family trees based on the estimated space-time ETAS model
"""
import numpy as np
# 1. For each pair of eartquakes i, j = 1,2,...,N (i<j), calculate the
# probability ρ_ij y φj from the final solution in algorithm 1.
p = []
f = []
N = 20
# 2. Set j=1
j = 1

while j < N + 1:
# 3. Generate a uniform random number Uj in [0,1].
    Uj = np.random.uniform(0,1,1)

    # 4. If Uj < φj, then the jth event is considered to be an immigrant
    #   (background event)

    immigrant = False

    if Uj < f[j]:
        immigrant = True
    else:
        immigrant = False

# 5. Otherwise, select the smallest I such that Uj< φj + sum_{i=1}^I ρ_ij
# Then the jth event is considered to be a descendant of the Ith event

# 6. If j = N, the terminate the algorithm; otherwise, set j = j+1 and
#   got o step 3.
    if j == N:
        pass
    else:
        j = j+1
        