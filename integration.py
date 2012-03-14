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

  g = exp(-1 * (d**2) / 2)

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
    if abs(ret[currenty] - u[currenty]) > 1 and flag:
      flag = False
      print "ERROR TOO BIG CHANGE"
  #we need to satisfy the conditions that f(-infty) = 1, f(infty) = 0, f(0) = 1/2
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
  if False and difference < 0:
    newret = ret.copy()
    for i, v in enumerate(ret):
      if i+difference < 0:
        newret[i] = 1
      else:
        newret[i] = ret[i+difference]
    ret = newret
  elif False:
    for i, v in enumerate(ret):
      if i+difference >= len(ret):
        ret[i] = 0
      else:
        ret[i] = ret[i+difference]
  return ret

def init_line():
  l.set_data([], [])
  return l,

def update_line(num, data, line):
  seconds = num * tresolution
  line.set_data(xvals, upoints[int(seconds/tresolution)])
  return line,

thetaresolution = .3
tresolution = 1

tstart = 0.0
tend = 5

thetastart = -pi/2
thetaend = pi/2
thetastart = -20
thetaend = 20

'''
if(len(sys.argv) > 4):
  tend = float(sys.argv[4])
if(len(sys.argv) > 3):
  xend = float(sys.argv[3])
if(len(sys.argv) > 2):
  tresolution = float(sys.argv[2])
if(len(sys.argv) > 1):
  xresolution = float(sys.argv[1])
'''

# the x-range of our algorithm
thetasize = int((thetaend-thetastart)/thetaresolution + 1)
thetapoints = arange(thetasize)
thetavals = list(float(x)*thetaresolution + thetastart for x in range(thetasize))


tsize = int((tend-tstart)/tresolution + 1)
tpoints = arange(tsize - 1)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))


upoints = zeros([tsize,thetasize])


# set the initial conditions,
# only one person infected
for i in range(thetasize):
  #upoints[0][i] = 1 if xvals[i] < 0 else 0
  power = map_fn(thetavals[i])
  cutoff = 500
  if power > cutoff:
    print "TOO BIG"
    power = cutoff
  elif power < -cutoff:
    print "TOO SMALL"
    power = -cutoff
  upoints[0][i] = 1.0 - 1.0/(1.0+exp(-power))

for t in tpoints[1:len(tpoints)]:
  print (t*tresolution)

  u = upoints[t-1]
  newu = integral(u, C)
  upoints[t] = newu

  plt.figure(t)
  plt.xlabel('X val')
  plt.title('Time graph of wave, iteration #' + str(t))

  x = thetavals
  y = upoints[t]

  plt.scatter(x, y)

'''
#We ignore any data outside of these bounds, for display purposes
xmin = -1000
xmax = 1000

min_index = 0
max_index = thetasize - 1

for i in range(thetasize):
  if map_fn(thetavals[i]) < xmax:
    max_index = i
  if map_fn(thetavals[i]) > xmin and min_index == 0:
    min_index = i

print "Min: %d, Max: %d" % (min_index, max_index)
print "Total Min: %d, Total Max: %d" % (0, thetasize-1)

#l, = plt.plot(xvals, upoints[0], 'r-')
plt.xlabel('theta')
plt.title('Time graph of wave, delta')
#line_ani = animation.FuncAnimation(fig, update_line, tsize, init_line, fargs=(upoints, l), interval=tresolution*500, blit=True)

  plt.figure(t)
  plt.xlabel('X val')
  plt.title('Time graph of wave, iteration #' + str(t))

  x = list(map_fn(v) for v in thetavals[min_index:max_index])
  y = upoints[t][min_index:max_index]

  plt.scatter(x, y)
plt.show()
plt.figure(2)
plt.xlim(xstart, xend)
plt.plot(xvals, upoints[1], 'r-')
plt.figure(3)
plt.plot(xvals, upoints[2], 'r-')
plt.figure(4)
plt.plot(xvals, upoints[3], 'r-')
plt.figure(5)
plt.plot(xvals, upoints[10], 'r-')
plt.show()
'''
