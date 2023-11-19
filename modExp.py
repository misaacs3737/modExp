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
  M = (n - r) // phi
  raisin = pow(a, phi, m) - 1
  sum = 0
  choose = 1
  Pc_exp = 1
  ell = 0
  for i in range(len(P)):
    et = E[i] // T[i]
    if ell<et:
      ell = et
  inverses = generateInversePairs(ell, m, P)
  for i in range(min(ell, M + 1)):
    sum = (sum + (choose * Pc_exp)) % m
    Pc_exp = (Pc_exp * raisin) % m
    choose = (((choose * (M - i)) % m) // inverses[i + 1][1] * inverses[i + 1][0]) % m
  R = pow(a, r, m)
  return (sum * R) % m
