from math import sin, exp
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


# This is the integral from 0 to N of:
# B(x-x')u(x',t) dx'
def integral(u):
  ret = u.copy()
  v = 1.85
  flag = True

  for currentx in range(xsize):
    #this total is the sum up to currentx of the whole formula
    total = 0.0

    for x in range(currentx):
      #this total is the sum of the inner integral
      total2 = 0.0
      xprimesize = xsize
      xprimeresolution = xresolution 
  
      for xprime in range(xprimesize):
      
        if xprime is 0 or xprime is (xprimesize - 1):
          coefficient2 = 1
        elif xprime % 2 is 0:
          coefficient2 = 2
        else:
          coefficient2 = 4

  
        total2 = total2 + coefficient2 * (B(x*xresolution-xprime*xprimeresolution) * u[xprime])
      total2 = total2 * xprimeresolution / 3

      if x is 0 or x is (xsize - 1):
        coefficient = 1
      elif xprime % 2 is 0:
        coefficient = 2
      else:
        coefficient = 4
      total += coefficient * (1 - u[x]) * total2

    total = total * xresolution / 3
    ret[currentx] = -total/v + 1
    if abs(ret[currentx] - u[currentx]) > 1 and flag:
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

xresolution = .5
tresolution = 1

tstart = 0.0
tend = 3

xstart = -30
xend = 30

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
xsize = int((xend-xstart)/xresolution + 1)
xpoints = arange(xsize)
xvals = list(float(x)*xresolution + xstart for x in range(xsize))


tsize = int((tend-tstart)/tresolution + 1)
tpoints = arange(tsize - 1)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))


upoints = zeros([tsize,xsize])

# set the initial conditions,
# only one person infected
for i in range(xsize):
  #upoints[0][i] = 1 if xvals[i] < 0 else 0
  upoints[0][i] = 1.0 - 1.0/(1.0+exp(-xvals[i]))


for t in tpoints:
  print (t*tresolution)

  first = 0
  last = xsize-1
  margin = .1
  for x in xpoints:
    if upoints[t][x] > 1-margin:
      first = x
    elif upoints[t][x] < margin and last == xsize - 1:
      last = x
  width = (last - first)
  reduction_factor = int(width/11)+1
  print "Reduction factor: %d" % reduction_factor

  u = upoints[t]
  newu = integral(u)

  if((t+1) != tsize):
    upoints[t+1] = newu

fig = plt.figure()
#l, = plt.plot(xvals, upoints[0], 'r-')
plt.xlim(xstart, xend)
plt.ylim(0, 1)
plt.xlabel('x')
plt.title('Time graph of wave')
#line_ani = animation.FuncAnimation(fig, update_line, tsize, init_line, fargs=(upoints, l), interval=tresolution*500, blit=True)
for t in tpoints:
  plt.plot(xvals, upoints[t])
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
