from numpy import arange, zeros
from scipy import polyfit 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

cutoff2 = 10
cutoff1 = 15

def left(data, threshold):
  l = 0
  for x in range(len(data)):
    if data[x] > (1-threshold):
      l = x
  return l

def right(data, threshold):
  r = len(data) - 1
  for x in reversed(range(len(data))):
    if data[x] < threshold:
      r = x
  return r

def plot_width(u):
  width = [0.0] * tsize
  for t in range(tsize):
    l = left(u[t], .05)
    r = right(u[t], .05)
    width[t] = r-l

  #x = tvals[int(cutoff1/tresolution):int(43/tresolution)]
  x = tvals[int(cutoff1/tresolution):]
  #y = width[int(cutoff1/tresolution):int(43/tresolution)]
  y = width[int(cutoff1/tresolution):]
  (m,b) = polyfit(x,y,1)
  print (b)
  print (m)
  plt.plot(x, y, 'r-')
  plt.xlim(cutoff1, tend)
  plt.show()

def plot_center(u):
  center = [0.0] * tsize
  for t in range(tsize):
    l = left(u[t], .05)
    r = right(u[t], .05)
    wave = u[t][l:r]
    mean = sum(wave) / (len(wave) if len(wave) > 0 else 1)
    center[t] = right(u[t], mean)

  #x = tvals[int(cutoff2/tresolution):int(45/tresolution)]
  x = tvals[int(cutoff2/tresolution):]
  #y = center[int(cutoff2/tresolution):int(45/tresolution)]
  y = center[int(cutoff2/tresolution):]
  (m,b) = polyfit(x,y,1)
  print (b)
  print (m)
  plt.plot(x, y, 'r-')
  plt.xlim(cutoff2, tend)
  plt.show()

#source = './datae20.500000-0.000000-1000.000000-0.500000-0.000000-50.000000.txt'
#values = source[8:(len(source)-4)].split('-')

'''
source_list = ['./data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.010000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.200000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.300000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.400000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.600000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c0.800000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.500000-c1.000000.txt']
'''

source_list = ['./data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.000000-c0.300000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.000000-c0.500000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.000000-c0.700000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.850000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.900000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.950000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.000000.txt']
                
'''
source_list = ['./data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r0.750000.txt',
               './data1.000000-0.000000-1000.000000-0.250000-0.000000-60.000000-r1.000000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.250000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.500000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r1.750000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r2.000000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r2.500000.txt',
               './data0.500000-0.000000-1000.000000-0.500000-0.000000-60.000000-r3.000000.txt']
'''
source = source_list[int(sys.argv[1])]
values = source[6:(len(source)-4)].split('-')

xresolution = float(values[0])
xstart = float(values[1])
xend = float(values[2])

xsize = int((xend-xstart)/xresolution + 1)
xpoints = arange(xsize)
xvals = list(float(x)*xresolution + xstart for x in range(xsize))


tresolution = float(values[3])
tstart = float(values[4])
tend = float(values[5])

tsize = int((tend-tstart)/tresolution + 1)
tpoints = arange(tsize)
tvals = list(float(t)*tresolution + tstart for t in range(tsize))

upoints = zeros([tsize,xsize])

t = 0
for line in open(source, 'r'):
  linevals = line.strip().split('\t')
  for x in range(xsize):
    upoints[t][x] = linevals[x]
  t += 1


def plot_animated():
  def init_line():
    l.set_data([], [])
    return l,

  def update_line(num, data, line):
    seconds = num * tresolution
    line.set_data(xvals, upoints[int(seconds/tresolution)])
    return line,

  fig = plt.figure()
  l, = plt.plot([], [], 'r-')
  plt.xlim(xstart, xend)
  plt.ylim(0, 1)
  plt.xlabel('x')
  plt.title('Time graph of wave')
  line_ani = animation.FuncAnimation(fig, update_line, tsize, init_line, fargs=(upoints, l), interval=tresolution*500, blit=True)
  #plt.show()

plot_center(upoints)
plot_width(upoints)
plot_animated()
