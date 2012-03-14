from math import sin, exp
from numpy import arange, zeros
from pylab import plot,xlabel,ylabel,show
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

c = float(sys.argv[1])
r = float(sys.argv[2])

def B(d):
  alpha = 8
  r_0 = 1
  p = (1 + (abs(d)*r_0)**(alpha))**(-1)
  e = c*exp(-r * abs(d))
  g = exp(-1 * (d**2) / 2)

  return e


# This is the integral from 0 to N of:
# B(x-x')u(x',t) dx'
def integral(u, reduction_factor):
  ret = u.copy()

  for x in range(xsize):
    # Evaluate the integral N different times
    # bounds are from a to b.
    total = 0.0
    xprimesize = xsize / reduction_factor
    xprimeresolution = xresolution * reduction_factor

    for xprime in range(xprimesize):
    
      if xprime is 0 or xprime is (xprimesize - 1):
        coefficient = 1
      elif xprime % 2 is 0:
        coefficient = 2
      else:
        coefficient = 4

      total = total + coefficient * (B(x*xresolution-xprime*xprimeresolution) * u[xprime*reduction_factor])
    ret[x] = total * xprimeresolution / 3
  return ret

def init_line():
  l.set_data([], [])
  return l,

def update_line(num, data, line):
  seconds = num * tresolution
  line.set_data(xvals, upoints[int(seconds/tresolution)])
  return line,

# This should be a function of u and t,
# but we already know u if we know x and t
def f(u, reduction_factor):
  return (1-u) * integral(u, reduction_factor)

xresolution = .5
tresolution = .5

tstart = 0.0
tend = 60

xstart = 0.0
xend = 1000

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
tpoints = arange(tsize)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))


upoints = zeros([tsize,xsize])

'''
read in data
'''
#source = './data1.000000-0.000000-200.000000-0.250000-0.000000-30.000000.txt'
#startdata = open(source, 'r').readlines()[40].strip().split('\t')

#for i in range(200):
  #upoints[0][i] = float(startdata[i])




# set the initial conditions,
# only one person infected
for i in range(xsize/9):
  upoints[0][i] = 1.0

data = open('./data%f-%f-%f-%f-%f-%f-r%f-c%f.txt' % (xresolution, xstart, xend, tresolution, tstart, tend, r, c), 'w')

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


  #for x in xpoints:
  u = upoints[t]
  k1 = tresolution*f(u, reduction_factor)
  k2 = tresolution*f(u+0.5*k1, reduction_factor)
  k3 = tresolution*f(u+0.5*k2, reduction_factor)
  k4 = tresolution*f(u+k3, reduction_factor)
  u += (k1+2*k2+2*k3+k4)/6

  for val in u:
    data.write('%e\t' % val)
  data.write('\n')

  if((t+1) != tsize):
    upoints[t+1] = u

data.close()

fig = plt.figure()
l, = plt.plot([], [], 'r-')
plt.xlim(xstart, xend)
plt.ylim(0, 1)
plt.xlabel('x')
plt.title('Time graph of wave')
line_ani = animation.FuncAnimation(fig, update_line, tsize, init_line, fargs=(upoints, l), interval=tresolution*500, blit=True)
plt.show()