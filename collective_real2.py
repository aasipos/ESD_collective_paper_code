#!/usr/bin/python
import sys
import numpy as np

def kopik(V1,V2,nu):
  tmp = V1**(1+nu)*V2**(1-nu)/(V1+V2)
  return tmp

if len(sys.argv) < 9:
  print "usage: %s N tmax v nu seed sample bins alpha" % ( sys.argv[0] )
  sys.exit(0)

N = int(sys.argv[1])
tmax = int(sys.argv[2])
v = float(sys.argv[3])
nu = float(sys.argv[4])
seed = int(sys.argv[5])
sample = int(sys.argv[6])
D = int(sys.argv[7])
alpha = float(sys.argv[8])
lr = np.zeros((N, sample),dtype=float)
np.random.seed(seed)
r = np.random.rand(N, sample) * 0.2 + 0.9
r[0] = alpha

i = np.zeros(sample,dtype=int)
j = np.zeros(sample,dtype=int)
dri = np.zeros(sample,dtype=float)
drj = np.zeros(sample,dtype=float)
for t in range(tmax):
  for tt in range(len(r)/2):
    [i, j] = np.random.choice(N,2,replace=False)
    dri = v * kopik(r[i],r[j],nu)
    drj = v * kopik(r[j],r[i],nu)
    r[i] -= dri
    r[j] -= drj
  for a in range(N):
    lr[a] = np.log(r[a])
  dist = np.zeros(D+1,dtype=int)
  ldist = np.zeros(D+1,dtype=int)
  Min = np.min(r)
  Max = np.max(r)
  lMin = np.min(lr)
  lMax = np.max(lr)
  for a in range(N):
    for b in range(sample):
      dist[int(D * (r[a][b] - Min) / (Max - Min))] += 1
      ldist[int(D * (lr[a][b] - lMin) / (lMax - lMin))] += 1
  dist[D - 1] += dist[D]
  ldist[D - 1] += ldist[D]
  for a in range(D):
    print a, (a + 0.5) * (Max - Min) / D + Min, dist[a], \
      (a + 0.5) * (lMax - lMin) / D + lMin, ldist[a], \
      Min, Max, lMin, lMax, np.mean(r), np.var(r), t, \
      np.mean(r[1:]), np.var(r[1:])
  print ""
  print ""
