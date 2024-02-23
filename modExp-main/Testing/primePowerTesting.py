import numpy as np
import math
from math import sqrt
import time
from timeit import default_timer as timer
import sys
import matplotlib.pyplot as plt
import sympy
import random

# compute the p-adic valuation of n
def vp(p, n):
  if (n % p == 0):
    return 1+vp(p, n//p)
  else:
    return 0

def inverse(i,m,L, p):
  inverses = L
  #finds the inverse of i mod m where L is the list of the inverses of 1 through i-1*
  #*when inversion is not posssible i.e. v = gcd(j, m) is not 1, then instead inverse of j/v is returned
  #v = math.gcd(i, m)
  v = 1
  for j in range(len(p)):
    v = v * p[j]**vp(p[j], i)
  if v == 1:
    inv = (inverses[m%i][0])*((m - m//(i))//(inverses[m%i][1]))%m
  else:
    inv = (inverses[(i)//v][0])
  inverses.append([inv, v])
  return inverses

def generateinverses(i, m, p):
  L = [[0,0], [1, 1]]
  for k in range(2, i+1):
    L = inverse(k,m,L, p)
  return L

def checkgoodvals(a,n,p,e,t):
  m = 1
  if len(p) != len(e) or len(e) != len(t):
    return "listsbad"
  for i in range(len(p)):
    if math.gcd(a,p[i])!= 1:
      return "a bad with: " + str(p[i])
    if t[i] > e[i]:
      return "t is big"
  return "good"

def binom_squareful(a,n,p,e, t):
  #m = product p[i]**e[i]
  T = 1
  phi = 1 #phi of T
  m = 1 #modulus
  p_min = min(p) #can be optimized
  for i in range(len(p)):
    t_i = p[i]**(t[i] - 1)
    phi = phi * t_i * (p[i] - 1)
    T = T * t_i * p[i]
    m = m * (p[i]**e[i])
  r = (n%phi)
  M = (n-r)//phi
  raisin = pow(a,phi,m)-1
  sum = 0
  choose = 1
  Pc_exp = 1
  et = []
  for i in range(len(p)):
    et.append(e[i]//t[i])
  ell = max(et)
  inverses = [0,1]
  v=1
  inverses = generateinverses(ell, m, p)
  for i in range(min(ell,M+1)):
    sum = (sum + (choose * Pc_exp)%m)%m
    Pc_exp = (Pc_exp * raisin)%m
    choose = (((choose * (M-i)) % m)//inverses[i+1][1] * inverses[i+1][0])%m
  R = pow(a,r,m)
  return (sum*R)%m


def silent_squareful_test(a,n,p,e,t,reps):
  start = timer()
  for i in range(reps):
    binom_ans = binom_squareful(a,n,p,e,t)
  end = timer()
  binom_time = end - start
  return [binom_ans, binom_time]

def silent_squareful_test_no_vp(a,n,p,e,t,reps):
  start = timer()
  for i in range(reps):
    binom_ans = binom_squareful(a,n,p,e,t)
  end = timer()
  binom_time = end - start
  return [binom_ans, binom_time]

def silent_powtest(a,n,p,k, reps):
  start = timer()
  pk=p**k
  for i in range(reps):
    powans = pow(a,n,pk)
  end = timer()
  powtime = end - start
  return [powans, powtime]

reps = 10
al = list(sympy.sieve.primerange(500000, 1000000))
p_list = [10000019] # initialize p
e_list = [int(np.log(p_list[0]))] # initialize e
tee = 1 # T = rad(n)
t_list = []
for i in range(len(p_list)):
  t_list.append(tee)


def run():
  m = p_list[0]**e_list[0]
  n = random.randrange(m//2, m) # Choose n in the interval [p^k/2, p^k]
  a = random.randrange(m//2, m) # Choose a in the interval [p^k/2, p^k]
  while a % p_list[0] == 0:
    a = random.randrange(m//2, m) # make sure gcd(a,p) = 1
  print(checkgoodvals(a,n,p_list, e_list, t_list))
  power = silent_powtest(a,n,m,1,reps)
  squareful = silent_squareful_test_no_vp(a,n,p_list, e_list,t_list, reps)
  powtime = power[1]
  squarefultime = squareful[1]
  dif = power[0] - squareful[0]
  return powtime/squarefultime

xpoints = []
ypoints = []

count = 0

for i in al:
  # iterate through prime list and run tests
  p_list[0] = i
  lg = int(np.log(p_list[0]))
  e_list[0] = random.randrange(int(lg-np.sqrt(lg)), int(lg+np.sqrt(lg))) # pick k randomly in the range [log p - sqrt(log p), log p + sqrt(log p)]
  ypoints.append(run())
  xpoints.append(count)
  count = count + 1

plt.plot(xpoints, ypoints)
plt.show()

