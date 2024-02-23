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
reps = 40
for i in range(10, 50):
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
  r.append(x[1]/x[0])
#print(y)
print(r)






#data analysis
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import linregress

r = [4.113797675688324, 3.324505068733082, 4.6763848396501455, 3.998389005786428, 5.452398824044536, 4.333511270297729, 6.543277082199238, 4.8953686229619775, 6.852990119007999, 5.779292082021989, 3.676639474303729, 8.331478312338357, 6.834877979523907, 6.610616233760443, 8.423175406701102, 7.036114958056762, 9.685008018217038, 7.289816788199973, 11.07865386496817, 7.969919876872934, 11.873773945485498, 8.553884299838284, 12.938796978561827, 9.48200348636744, 12.965260448551222, 9.694313901143754, 15.13615272239067, 9.82457431439049, 14.83834697761794, 15.10710571236946, 16.43075935255352, 11.07514699874151, 17.222834028176177, 12.26108472605459, 17.999161798598756, 18.05248075718227, 18.515366136071833, 13.321269473757505, 20.313727758466012, 13.869658004187086]

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



def sqrt_fit(x, c):
    return c * np.sqrt(x)

# Perform curve fitting for odd and even lists
params_even, _ = curve_fit(sqrt_fit, logmeven, reven)
params_odd, _ = curve_fit(sqrt_fit, logmodd, rodd)

# Generate fitted curves
x_values = np.linspace(min(logm), max(logm), 1000)
fitted_even = sqrt_fit(x_values, *params_even)
fitted_odd = sqrt_fit(x_values, *params_odd)

# Plot the data and fitted curves
plt.scatter(logmodd, rodd, c="orange", label="Odd n", s=50)
plt.scatter(logmeven, reven, c="blue", label="Even n", s=50)
plt.plot(x_values, fitted_even, color="blue", linestyle="--", label="Even fit: c=%.3f" % params_even[0])
plt.plot(x_values, fitted_odd, color="orange", linestyle="--", label="Odd fit: c=%.3f" % params_odd[0])
plt.xlabel("log m")
plt.ylabel("Ratio")
plt.legend()
plt.show()



# Calculate residuals for even and odd fits
residuals_even = reven - sqrt_fit(logmeven, *params_even)
residuals_odd = rodd - sqrt_fit(logmodd, *params_odd)

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
