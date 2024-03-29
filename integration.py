from math import sin, exp, tan, atanh, pi, cos, floor
from numpy import arange, zeros, convolve
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
  return 1/(cos(theta)**2)

def gaussian1(d):
  return exp(-1 * ((.125*d)**2) / 2)

def gaussian2(d):
  return exp(-1 * ((.25*d)**2) / 2)

def gaussian3(d):
  return exp(-1 * ((.5*d)**2) / 2)

def gaussian4(d):
  return exp(-1 * ((1*d)**2) / 2)

def gaussian5(d):
  return exp(-1 * ((2*d)**2) / 2)

def delta1(d):
  if abs(d) < .00001:
    return 3 / thetaresolution
  else:
    return 0


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

  #precalculate convolution:
  conv = [0.0] * thetasize
  for y in range(thetasize):
    conv[y] = 0.0
    for yprime in range(thetasize):
      coeff2 = compute_coeff(yprime, thetasize)
      conv[y] += coeff2 * (beta(map_fn(thetavals[y])-map_fn(thetavals[yprime])) * u[yprime]) * d_map_fn(thetavals[yprime])
    conv[y] = conv[y] * thetaresolution / 3

  for currenty in range(thetasize):
    #this total is the sum up to currentx of the whole formula
    total = 0.0

    #beta_vals = 
    #print convolve(beta_vals, u, 'same')
    for y in range(currenty):
      #this total is the sum of the inner integral
      coefficient = compute_coeff(y, thetasize)
      total += coefficient * (1 - u[y]) * conv[y] * d_map_fn(thetavals[y])

    total = total * thetaresolution / 3
    ret[currenty] = -total

  #we need to satisfy the conditions that f(-infty) = 1, f(infty) = 0, f(0) = 1/2
  return rescale(ret)

def rescale(ret):
  #first scale ret
  v = abs(ret[len(ret)-1])
  print "Velocity: " + str(v)
  if v <= 0:
    return ret
  ret = ret/v
  #then shift ret upward
  ret = ret + 1
  #then shift sideways
  middle_index = 0
  expected_middle_index = floor(len(ret)/2)
  for i, v in enumerate(ret):
    if v < .5:
      middle_index = i
      break
  difference = middle_index - expected_middle_index
  print "Difference: %d" % difference
  if difference < 0:
    newret = ret.copy()
    for i, v in enumerate(ret):
      if i+difference < 0:
        newret[i] = 1
      else:
        newret[i] = ret[i+difference]
    ret = newret
  else:
    for i, v in enumerate(ret):
      if i+difference >= len(ret):
        ret[i] = 0
      else:
        ret[i] = ret[i+difference]
  return ret

thetaresolution = .2
tresolution = 1

tstart = 0.0
tend = 8

thetastart = -40
thetaend = 40


# the x-range of our algorithm
thetasize = int((thetaend-thetastart)/thetaresolution + 1)
thetapoints = arange(thetasize)
thetavals = list(float(x)*thetaresolution + thetastart for x in range(thetasize))

#We ignore any data outside of these bounds, for display purposes
xmin = -100
xmax = 100

min_index = 0
max_index = thetasize - 1

for i in range(thetasize):
  if map_fn(thetavals[i]) < xmin and min_index == 0:
    max_index = i
  if map_fn(thetavals[i]) > xmax:
    min_index = i

print "Min: %d, Max: %d" % (min_index, max_index)
print "Total Min: %d, Total Max: %d" % (0, thetasize-1)

thetavals = thetavals[min_index:max_index]
thetapoints = thetapoints[min_index:max_index]
thetasize = len(thetavals)

tsize = int((tend-tstart)/tresolution + 1)
tpoints = arange(tsize - 1)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))

upoints1 = zeros([tsize,thetasize])
upoints2 = zeros([tsize,thetasize])
upoints3 = zeros([tsize,thetasize])
upoints4 = zeros([tsize,thetasize])
upoints5 = zeros([tsize,thetasize])


# set the initial conditions,
# only one person infected
for i in range(thetasize):
  power = map_fn(thetavals[i])

  upoints1[0][i] = 1.0 - 1.0/(1.0+exp(-power/5))
  upoints2[0][i] = 1.0 - 1.0/(1.0+exp(-power/4))
  upoints3[0][i] = 1.0 - 1.0/(1.0+exp(-power/3))
  upoints4[0][i] = 1.0 - 1.0/(1.0+exp(-power/2))
  upoints5[0][i] = 1.0 - 1.0/(1.0+exp(-power/1))

for t in tpoints[1:len(tpoints)]:
  print (t*tresolution)

  u = upoints1[t-1]
  newu = integral(u, gaussian1)
  upoints1[t] = newu

  u = upoints2[t-1]
  #newu = integral(u, gaussian2)
  upoints2[t] = newu

  u = upoints3[t-1]
  newu = integral(u, gaussian1)
  upoints3[t] = newu

  u = upoints4[t-1]
  #newu = integral(u, gaussian4)
  upoints4[t] = newu
  
  u = upoints5[t-1]
  newu = integral(u, gaussian1)
  upoints5[t] = newu

for t in tpoints:
  plt.figure(t)
  plt.xlabel('X val')
  plt.title('Time graph of wave, iteration #' + str(t))

  x1 = list(map_fn(v) for v in thetavals)
  x3 = list(map_fn(v)*2 for v in thetavals)
  x5 = list(map_fn(v)*4 for v in thetavals)

  plt.scatter(x1, upoints1[t], color = 'b')
  #plt.scatter(x, upoints2[t], color = 'r')
  plt.scatter(x1, upoints3[t], color = 'g')
  #plt.scatter(x, upoints4[t], color = 'k')
  plt.scatter(x1, upoints5[t], color = 'y')

plt.figure(len(tpoints))
plt.xlabel('Distance')
plt.title('Beta functions')

plt.scatter(x1, list(gaussian1(map_fn(x)) for x in thetavals), color = 'b')
#plt.scatter(x, list(gaussian2(map_fn(x)) for x in thetavals), color = 'r')
plt.scatter(x1, list(gaussian1(map_fn(x)) for x in thetavals), color = 'g')
#plt.scatter(x, list(gaussian4(map_fn(x)) for x in thetavals), color = 'k')
plt.scatter(x1, list(gaussian1(map_fn(x)) for x in thetavals), color = 'y')
plt.show()
