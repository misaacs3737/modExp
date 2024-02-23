import time
import matplotlib.pyplot as plt
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


#testing
P = [101]
E = [200]
T = [13]
k = len(P)
m = 1
for i in range(k):
        m *= P[i]**E[i]
a = 13
n = m//3
times = []
Tlist = []
for i in range(1, 50, 1):
    T = [i]
    s = time.time()
    for j in range(30):
        x = ourModExp(a, n, P, E, T)
    e = time.time()
    ourtime = e - s
    times.append(ourtime/30)
    Tlist.append(i)
plt.plot(Tlist, times)
# time goes as at + b/t, need to solve for a,b by plugging in say t = t1, t = t2
t1 = 30
t2 = 8
T1 = times[t1]
T2 = times[t2]
d = 1/(t1/t2 - t2/t1)
a = d * (T1/t2 - T2/t1)
b = d * (-1*T1*t2 + t1*T2)
pred = []
for i in range(1, 50, 1):
   pred.append(a*i + b/i)
plt.plot(Tlist, pred)
plt.xlabel("t value")
plt.ylabel("time (s)")
plt.show()
