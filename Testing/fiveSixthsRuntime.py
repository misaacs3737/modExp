#data collection

import time

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

def times(reps, a , n , P, E, T): #returns [ours, builtin] or [-1, -1] if incorrect computation
    start = time.time()
    for i in range(reps):
        ournum = ourModExp(a, n, P, E, T)
    end = time.time()
    ourtime = end - start
    start = time.time()
    for i in range(reps):
        m=1
        for j in range(len(P)):
          m *= P[j]**E[j]
        builtinnum = pow(a, n, m)
    end = time.time()
    builtintime = end - start
    if builtinnum != ournum:
      return [-1, -1]
    return [ourtime, builtintime]
#51 elements - primes[i] ~= 10^i
y = []
r = []

primes = [2, 11, 101,  1009, 10007, 100003,1000003,10000019,100000007,1000000007,10000000019,100000000003, 1000000000039,10000000000037, 100000000000031, 1000000000000037, 10000000000000061, 100000000000000003, 1000000000000000003, 10000000000000000051, 100000000000000000039, 1000000000000000000117,10000000000000000000009, 100000000000000000000117,1000000000000000000000007,10000000000000000000000013, 100000000000000000000000067, 1000000000000000000000000103, 10000000000000000000000000331, 100000000000000000000000000319, 1000000000000000000000000000057, 10000000000000000000000000000033, 100000000000000000000000000000049,1000000000000000000000000000000061, 10000000000000000000000000000000193,100000000000000000000000000000000069,1000000000000000000000000000000000067,10000000000000000000000000000000000043, 100000000000000000000000000000000000133,1000000000000000000000000000000000000003, 10000000000000000000000000000000000000121, 100000000000000000000000000000000000000109, 1000000000000000000000000000000000000000063, 10000000000000000000000000000000000000000057, 100000000000000000000000000000000000000000031, 1000000000000000000000000000000000000000000009, 10000000000000000000000000000000000000000000121, 100000000000000000000000000000000000000000000033, 1000000000000000000000000000000000000000000000193, 10000000000000000000000000000000000000000000000009, 100000000000000000000000000000000000000000000000151]
reps = 50
delvalues = 40
for i in range(10, 10+delvalues):
  #p ~= 10^(i)
  P = [primes[i]]
  T = [1]
  E = [i] # log_10 p  
  a = 13424
  n = (P[0]**E[0] * 2) // 3
  x = times(reps, a, n, P, E, T)
  if x==[-1, -1]:
    print(i)
  y.append(x)
  
  r.append(x[1]**(-5/6)*x[0])
#print(y)
print(r)


#data analysis
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import linregress

r = [0.12237074092416536, 0.17060494458136372, 0.0894578954909992, 0.16164029979290276, 0.11297551375789891, 0.14787121196635175, 0.1027909908120211, 0.14323825153822964, 0.10649125709705247, 0.1449775745335321, 0.1457176407238139, 0.10488820799160042, 0.13102338854032916, 0.14408886225256093, 0.0985682573032332, 0.14896575302600337, 0.10360528354234926, 0.1287729748650683, 0.11400562515791358, 0.1495235343218217, 0.0952630928925638, 0.13691667760103332, 0.09569071139347235, 0.14056599725153007, 0.10166491250801409, 0.1351806309007308, 0.09169531628399573, 0.13958396781488824, 0.09933362168205649, 0.10289179502955183, 0.09161266366360521, 0.13936146025443172, 0.09648680664803143, 0.13517226098217308, 0.09742755018587874, 0.09656556101780737, 0.09459459724073363, 0.13669969547720795, 0.09442257253684964, 0.13685554770686711]
logm = []
for i in range(10, 50):
    logm.append((i)**2)
reven = []
rodd = []
logmeven = []
logmodd = []

for i in range(0, 40):
    if i%2!=0:
        logmodd.append(logm[i])
        rodd.append(r[i])
    else:
        logmeven.append(logm[i])
        reven.append(r[i])
def average(L):
    s = 0
    for i in range(len(L)):
        s+= L[i]
    return s/len(L)

oddavg = average(rodd)
evenavg = average(reven)
print(oddavg/evenavg)
#fill in this code



def line_fit(x, c):
    return c * np.sqrt(x)/np.sqrt(x)

# Perform curve fitting for odd and even lists
params_even, _ = curve_fit(line_fit, logmeven, reven)
params_odd, _ = curve_fit(line_fit, logmodd, rodd)

# Generate fitted curves
x_values = np.linspace(min(logm), max(logm), 1000)
fitted_even = line_fit(x_values, *params_even)
fitted_odd = line_fit(x_values, *params_odd)

# Plot the data and fitted curves
plt.scatter(logmodd, rodd, c="orange", label="Odd i", s=50)
plt.scatter(logmeven, reven, c="blue", label="Even i", s=50)
plt.plot(x_values, fitted_even, color="blue", linestyle="--", label="Even fit: c=%.3f" % params_even[0])
plt.plot(x_values, fitted_odd, color="orange", linestyle="--", label="Odd fit: c=%.3f" % params_odd[0])
plt.xlabel("log m")
plt.ylabel("T₀T₁^(-5/6)")
plt.legend()
plt.show()



# Calculate residuals for even and odd fits
residuals_even = reven - line_fit(logmeven, *params_even)
residuals_odd = rodd - line_fit(logmodd, *params_odd)

# Calculate R-squared for even fit
ss_total_even = np.sum((reven - np.mean(reven))**2)
ss_res_even = np.sum(residuals_even**2)
r_squared_even = 1 - (ss_res_even / ss_total_even)

# Calculate R-squared for odd fit
ss_total_odd = np.sum((rodd - np.mean(rodd))**2)
ss_res_odd = np.sum(residuals_odd**2)
r_squared_odd = 1 - (ss_res_odd / ss_total_odd)

print("R-squared for even fit:", r_squared_even)
print("R-squared for odd fit:", r_squared_odd)
