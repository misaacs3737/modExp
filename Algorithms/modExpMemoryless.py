# returns the exponent of prime p in the prime factorization of n
def nu(p, n):
  if (n % p == 0):
    return 1+nu(p, n // p)
  else:
    return 0
# returns a^n modulo m = product P[i]**E[i]
def ourModExpMemoryless(a , n , P, E, T):
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
  for i in range(min(ell, q + 1)):
    sum = (sum + (choose * cExp)) % m
    cExp = (cExp * c) % m
    v = 1
    for j in range(len(P)):
      v *= P[j]**nu(P[j], i+1)
    u = pow((i+1)//v, -1, m)
    choose = (((choose * (q - i)) % m) // v * u) % m
  ar = pow(a, r, m)
  return (sum * ar) % m
