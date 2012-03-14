from math import sin, exp, tan, pi
from numpy import arange, zeros
from pylab import plot,xlabel,ylabel,show
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

def B(d):
  alpha = 8
  r_0 = 1
  c = 1
  r = 1
  p = (1 + (abs(d)*r_0)**(alpha))**(-1)
  e = c*exp(-r * abs(d))
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
def integral(u):
  ret = u.copy()
  v = 1.85
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
        coefficient2 = compute_coeff(yprime, yprimesize)
        total2 = total2 + coefficient2 * (B(tan(y*thetaresolution)-tan(yprime*yprimeresolution)) * u[yprime])
      total2 = total2 * yprimeresolution / 3

      coefficient = compute_coeff(y, thetasize)
      total += coefficient * (1 - u[y]) * total2

    total = total * thetaresolution / 3
    ret[currenty] = -total/v + 1
    if abs(ret[currenty] - u[currenty]) > 1 and flag:
      flag = False
      print "ERROR TOO BIG CHANGE"
  return ret

def init_line():
  l.set_data([], [])
  return l,

def update_line(num, data, line):
  seconds = num * tresolution
  line.set_data(xvals, upoints[int(seconds/tresolution)])
  return line,

thetaresolution = .05
tresolution = 1

tstart = 0.0
tend = 2

thetastart = -pi/2
thetaend = pi/2

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
  power = tan(thetavals[i])
  cutoff = 500
  if power > cutoff:
    power = cutoff
  elif power < -cutoff:
    power = -cutoff
  print power
  upoints[0][i] = 1.0 - 1.0/(1.0+exp(-power))


for t in tpoints:
  print (t*tresolution)

  u = upoints[t]
  newu = integral(u)

  if((t+1) != tsize):
    upoints[t+1] = newu

fig = plt.figure()
#l, = plt.plot(xvals, upoints[0], 'r-')
plt.xlim(thetastart, thetaend)
plt.ylim(0, 1)
plt.xlabel('theta')
plt.title('Time graph of wave')
#line_ani = animation.FuncAnimation(fig, update_line, tsize, init_line, fargs=(upoints, l), interval=tresolution*500, blit=True)
for t in tpoints:
  plt.plot(thetavals, upoints[t])
  plt.show()

'''
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
