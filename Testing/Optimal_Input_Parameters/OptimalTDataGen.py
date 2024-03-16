import time
import matplotlib.pyplot as plt
import numpy as np
import sympy
from math import log, sqrt

# returns the exponent of prime p in the prime factorization of n
def nu(p, n):
  if (n % p == 0):
    return 1+nu(p, n // p)
  else:
    return 0

# returns the inverse pair of i, given the inverse pairs of 0 through i-1, inclusive
def nextInversePair(i,m,L, P):
  v = 1
  for j in range(len(P)):
    v *= P[j]**nu(P[j], i)
  if v == 1:
    u = (L[m % i][0]) * ((m - m // i) // (L[m % i][1])) % m
  else:
    u = (L[i // v][0])
  return [u, v]

# returns the inverse pairs of 0 through i, inclusive
def generateInversePairs(i, m, p):
  L = [[0, 0], [1, 1]]
  for j in range(2, i+1):
    L.append(nextInversePair(j, m, L,  p))
  return L

# returns a^n modulo m = product P[i]**E[i]
def ourModExp(a , n , P, E, T):
  t = 1
  phi = 1
  m = 1
  for i in range(len(P)):
    temp = P[i]**(T[i] - 1)
    phi *= temp * (P[i] - 1)
    t *= temp * P[i]
    m *= P[i]**E[i]
  r = n % phi
  q = (n - r) // phi
  c = pow(a, phi, m) - 1
  sum = 0
  choose = 1
  cExp = 1
  ell = 0
  for i in range(len(P)):
    et = E[i] // T[i]
    if ell<et:
      ell = et
  inverses = generateInversePairs(ell, m, P)
  for i in range(min(ell, q + 1)):
    sum = (sum + (choose * cExp)) % m
    cExp = (cExp * c) % m
    choose = (((choose * (q - i)) % m) // inverses[i + 1][1] * inverses[i + 1][0]) % m
  ar = pow(a, r, m)
  return (sum * ar) % m

# computes optimal t input
def compute_t(a,n,p,k):
  P = [p]
  E = [k]
  optimal_t = 1
  curr_time = 100000
  for i in range(1, 50, 1):
    T = [i]
    s = time.time()
    for j in range(3):
      x = ourModExp(a,n,P,E,T)
    e = time.time()
    if (e-s<curr_time):
      curr_time = e - s
      optimal_t = i
  return optimal_t

# primes between p_1500 and p_3500
primes = list(sympy.primerange(1500, 3500))

# Generate data and compute t
data = []
count = 1
for p in primes:
    a = np.random.randint(1, p//4) # choosing a randomly between 1 and p/4
    k = 30 + np.random.randint(log(p) - sqrt(log(p)), log(p) + sqrt(log(p))) # choosing k randomly in a logarithmic interval with a buff of 30
    n = np.random.randint(3,10)
    n = (p ** k)//n # choosing n randomly between p^k/10 and p^k/3 due to randint overflows
    t = compute_t(a, n, p, k)
    data.append([a, p, k, t]) # omitting n from the dataset due to size reasons
    print(count) # to track generation at runtime
    count = count + 1

# Convert data to numpy array
data = np.array(data)

# Save data to a file
start = time.time()
np.save('optimal_t_data.npy', data)
end = time.time()
print("Finished in "  + str(end-start))

