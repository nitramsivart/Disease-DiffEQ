from math import sin, exp, tan, atanh, pi, cos, floor
from numpy import arange, zeros, convolve
from pylab import plot,xlabel,ylabel,show
from multiprocessing import Pool
from scipy import linalg, array
import matplotlib.pyplot as plt
import matplotlib.animation as animation

thetaresolution = .4
tresolution = .2

tstart = 0.0
tend = 20

thetastart = -30
thetaend = 30

sign = -1

velocity = 40
start_val = .5

def initial_condition(power):
  if abs(power) < .5:
    return exp(-1./(1-power**2))/20
  else:
    return 0.
  return 1.0 - 1.0/(1.0+exp(-power))
  '''
  upoints1[0][i] = start_val
  upoints1[0][i] = 1.0 - float(i)/(thetasize - 1)
  upoints1[0][i] = 1.0 - 1.0/(1.0+exp(-power))
  '''


def map_fn(theta):
  return theta
  return tan(theta)

def d_map_fn(theta):
  return 1
  return 1/(cos(theta)**2)

def gaussian1(d):
  return exp(-1 * ((d)**2) / 2)

def exp1(d):
  return exp(-1 * ( d / 2))

def exp2(d):
  return exp(-1 * (.5 * d / 2))

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

def power1(d):
  r_0 = 1
  alpha = 1
  return (1 + (abs(d)*r_0)**(alpha))**(-1)
  return d**(-1)

def compute_coeff(index, size):
  if index is 0 or index is (size - 1):
    return 1
  elif index % 2 is 0:
    return 2
  else:
    return 4

def compute_deriv(f, xindex, xwidth):
  if xindex == (len(f) - 1):
    return (f[xindex] - f[xindex - 1])/xwidth
  elif xindex == 0:
    return (f[xindex + 1] - f[xindex])/xwidth
  else:
    return (f[xindex + 1] - f[xindex - 1]) / (2 * xwidth)

# This is the integral from 0 to N of:
# B(x-x')u(x',t) dx'
def integral(u, beta, vel):
  change = u.copy()
  flag = True

  #precalculate convolution:
  conv = [0.0] * thetasize
  for y in range(thetasize):
    conv[y] = 0.0
    for yprime in range(thetasize):
      coeff2 = compute_coeff(yprime, thetasize)
      conv[y] += coeff2 * (beta(map_fn(thetavals[y])-map_fn(thetavals[yprime])) * u[yprime]) * d_map_fn(thetavals[yprime])
    conv[y] = conv[y] * thetaresolution / 3
  
  #calculate new values using implicit method
  h = tresolution
  a = thetaresolution

  A = array([[h/(2*a)]*(thetasize), [sign * 1.]*(thetasize), [-h/(2*a)]*(thetasize)])
  v = (array([h*(-1/float(vel))*(1-u[x])*conv[x] + sign * u[x] for x in range(thetasize)])).T
  v[0] = v[0] + h/(2*a) * u[0]
  v[thetasize-1] = v[thetasize-1] - h/(2*a) * u[thetasize-1]

  sol = linalg.solve_banded((1,1), A, v)
  return sol

def do_runge_kutta(u, beta, v):
  return integral(u, beta, v)

  k1 = integral(u, beta, v)
  k2 = integral(u+.5*k1, beta, v)
  k3 = integral(u+.5*k2, beta, v)
  k4 = integral(u+k3, beta, v)

  ret = u.copy()
  ret += (k1+2*k2+2*k3+k4)/6

  return ret
  return rescale(ret)

def rescale(ret):
  '''
  #first scale ret
  v = abs(ret[len(ret)-1])
  print "Velocity: " + str(v)
  if v <= 0:
    return ret
  ret = ret/v
  #then shift ret upward
  ret = ret + 1
  '''
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


# the x-range of our algorithm
thetasize = int((thetaend-thetastart)/thetaresolution + 1)
thetapoints = arange(thetasize)
thetavals = list(float(x)*thetaresolution + thetastart for x in range(thetasize))

#We ignore any data outside of these bounds, for display purposes
xmin = -200
xmax = 200

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

def init_line1():
  l1.set_data([], [])
  return l1,

def update_line1(num, data, line):
  line.set_data(thetavals, upoints1[num])
  return line,

def init_line2():
  l2.set_data([], [])
  return l2,

def update_line2(num, data, line):
  line.set_data(thetavals, upoints2[num])
  return line,

# set the initial conditions,
# only one person infected
for i in range(thetasize):
  power = map_fn(thetavals[i])
  upoints1[0][i] = initial_condition(power)


for t in tpoints[1:len(tpoints)]:
  print (t*tresolution)

  u = upoints1[t-1]
  newu = do_runge_kutta(u, gaussian1, velocity)
  upoints1[t] = newu

  u = upoints2[t-1]
  #newu = do_runge_kutta(u, gaussian1, 140)
  upoints2[t] = newu

  u = upoints3[t-1]
  #newu = integral(u, gaussian3)
  upoints3[t] = newu

  u = upoints4[t-1]
  #newu = integral(u, gaussian4)
  upoints4[t] = newu
  
  u = upoints5[t-1]
  #newu = integral(u, gaussian5)
  upoints5[t] = newu

fig1 = plt.figure(0)
plt.xlabel('X val')
plt.axis([thetastart, thetaend, 0, 1])
plt.title('Time graph of wave (power)' + str(velocity) + ', ' + str(start_val))
l1, = plt.plot(thetavals, upoints1[0], 'r')
line_ani = animation.FuncAnimation(fig1, update_line1, tsize, init_line1, fargs=(upoints1, l1), interval=tresolution*100, blit=True)

#fig2 = plt.figure(1)
#plt.xlabel('X val')
#plt.title('Time graph of wave' + str(thetaresolution) + ',' + str(tresolution) + ',' + str(tend))
#l2, = plt.plot(thetavals, upoints2[0], 'b')
#line_ani2 = animation.FuncAnimation(fig2, update_line2, tsize, init_line2, fargs=(upoints2, l2), interval=tresolution*50, blit=True)

plt.figure(len(tpoints))
plt.xlabel('Distance')
plt.title('Beta functions')

plt.scatter(thetavals, list(gaussian1(abs(map_fn(x))) for x in thetavals), color = 'b')
#plt.scatter(thetavals, list(power1(abs(map_fn(x))) for x in thetavals), color = 'r')

plt.show()
