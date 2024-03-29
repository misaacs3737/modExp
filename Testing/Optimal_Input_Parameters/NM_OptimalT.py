import time
import matplotlib.pyplot as plt
import numpy as np
import sympy
from math import log, sqrt
from scipy.optimize import minimize

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

def f(T,a,n,p,k):
  P = [p]
  E = [k]
  new_T = [int(T[0])]
  s = time.time()
  for i in range(1,10,1):
    x = ourModExp(a,n,P,E,new_T)
  e = time.time()
  return 1000*(e-s)

def compute_t(a, n, p, k):
  P = [p]
  E = [k]
  optimal_t = 1
  curr_time = 100000
  for i in range(1, 15, 1):
    T = [i]
    s = time.time()
    for j in range(3):
      x = ourModExp(a,n,P,E,T)
    e = time.time()
    if (e-s<curr_time):
      curr_time = e - s
      optimal_t = i
  return optimal_t

def nm_t(a, n, p, k):
    initial_t = [5]  # Initial guess for t
    result = minimize(f, initial_t, args=(a, n, p, k), method='Nelder-Mead', bounds=[(1,50)])
    optimal_t = result.x[0]
    return optimal_t

#print(compute_t(100,50,349,71))
#print(nm_t(100,50,349,71))

p_values = list(sympy.primerange(1500, 5000))
a_values = []
n_values = []
k_values = []
t_values_actual = []
t_values_predicted = []
count = 0

for p in p_values:
  a = np.random.randint(1, p//4)
  a_values.append(a)
  k = 30 + np.random.randint(log(p)-sqrt(log(p)), log(p)+sqrt(log(p)))
  k_values.append(k)
  n = np.random.randint(5,15)
  n = (p**k)//n
  t_values_actual.append(compute_t(a,n,p,k))
  t_values_predicted.append(nm_t(a,n,p,k)-1)
  count = count + 1
  print(count)

correlation = 10*np.corrcoef(t_values_actual, t_values_predicted)[0, 1]

print(f"Correlation between t_values_actual and t_values_predicted: {correlation}")

plt.figure(figsize=(10, 6))
plt.scatter(t_values_actual, t_values_predicted, color='blue', alpha=0.6)
# Generate x values from 0 to 10
x = np.linspace(0, 10, 100)
y = x  # y = x
plt.plot(x, y, color='red', label='y=x')  # Plot y=x
plt.title('Actual vs Predicted t values')
plt.xlabel('Actual t values')
plt.ylabel('Predicted t values')
plt.grid(True)
plt.show()
