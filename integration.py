from math import sin, exp, tan, atanh, pi, cos, floor
from numpy import arange, zeros
from pylab import plot,xlabel,ylabel,show
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

def map_fn(theta):
  return theta
  return tan(theta)

def d_map_fn(theta):
  return 1
  return theta
  return 1/(cos(theta)**2)

def gaussian1(d):
  return exp(-1 * (d**2) / 2)

def gaussian2(d):
  if abs(d) < .00001:
    return 3 / thetaresolution
  else:
    return 0

def gaussian3(d):
  return exp(-1 * (d**2) / 18)

def C(d):
  alpha = 8
  r_0 = 1
  c = 1
  r = 1
  p = (1 + (abs(d)*r_0)**(alpha))**(-1)
  e = c*exp(-r * abs(d))

  if abs(d) < .00001:
    return 1
  else:
    return 0


  return g

def compute_coeff(index, size):
  if index is 0 or index is (size - 1):
    return 1
  elif index % 2 is 0:
    return 2
  else:
    return 4

# This is the integral from 0 to N of:
# B(x-x')u(x',t) dx'
def integral(u, beta):
  ret = u.copy()
  flag = True

  for currenty in range(thetasize):
    #this total is the sum up to currentx of the whole formula
    total = 0.0

    for y in range(currenty):
      #this total is the sum of the inner integral
      total2 = 0.0
      yprimesize = thetasize
      yprimeresolution = thetaresolution
  
      for yprime in range(yprimesize):
        coeff2 = compute_coeff(yprime, yprimesize)
        total2 += coeff2 * (beta(map_fn(thetavals[y])-map_fn(thetavals[yprime])) * u[yprime]) * d_map_fn(thetavals[yprime])
      total2 = total2 * yprimeresolution / 3

      coefficient = compute_coeff(y, thetasize)
      total += coefficient * (1 - u[y]) * total2 * d_map_fn(thetavals[y])

    total = total * thetaresolution / 3
    ret[currenty] = -total
  #we need to satisfy the conditions that f(-infty) = 1, f(infty) = 0, f(0) = 1/2
  #first scale ret
  v = abs(ret[len(ret)-1])
  print "Velocity: " + str(v)
  if v <= 0:
    return ret
  ret = ret/v
  #then shift ret upward
  ret = ret + 1
  return ret


thetaresolution = .4
tresolution = 1

tstart = 0.0
tend = 6

thetastart = -10
thetaend = 10


# the x-range of our algorithm
thetasize = int((thetaend-thetastart)/thetaresolution + 1)
thetapoints = arange(thetasize)
thetavals = list(float(x)*thetaresolution + thetastart for x in range(thetasize))


tsize = int((tend-tstart)/tresolution + 1)
tpoints = arange(tsize - 1)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))


upoints1 = zeros([tsize,thetasize])
upoints2 = zeros([tsize,thetasize])
upoints3 = zeros([tsize,thetasize])


# set the initial conditions,
# only one person infected
for i in range(thetasize):
  upoints1[0][i] = 1.0 - 1.0/(1.0+exp(-thetavals[i]))
  upoints2[0][i] = 1.0 - 1.0/(1.0+exp(-thetavals[i]))
  upoints3[0][i] = 1.0 - 1.0/(1.0+exp(-thetavals[i]))

plt.figure(0)
plt.xlabel('X val')
plt.title('Time graph of wave, iteration #' + str(0))

plt.scatter(thetavals, upoints1[0])
plt.scatter(thetavals, upoints2[0], color = 'r')
#plt.scatter(thetavals, upoints3[0], color = 'g')

for t in tpoints[1:len(tpoints)]:
  print (t*tresolution)

  u = upoints1[t-1]
  newu = integral(u, gaussian1)
  upoints1[t] = newu

  u = upoints2[t-1]
  newu = integral(u, gaussian2)
  upoints2[t] = newu

  u = upoints3[t-1]
  #newu = integral(u, gaussian3)
  upoints3[t] = newu

  plt.figure(t)
  plt.xlabel('X val')
  plt.title('Time graph of wave, iteration #' + str(t))

  plt.scatter(thetavals, upoints1[t])
  plt.scatter(thetavals, upoints2[t], color = 'r')
  #plt.scatter(thetavals, upoints3[t], color = 'g')

plt.figure(len(tpoints))
plt.xlabel('Distance')
plt.title('Beta functions')
plt.scatter(thetavals, list(gaussian1(x) for x in thetavals))
plt.scatter(thetavals, list(gaussian2(x) for x in thetavals), color = 'r')
#plt.scatter(thetavals, gaussian3(thetavals), color = 'g')
plt.show()
