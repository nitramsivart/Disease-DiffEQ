from pylab import *

velocity = [19.2, 7.177, 5.184, 3.944, 2.923, 1.77, 1.27]
width = [38.33, 22.5, 18.7, 16.0, 14.16, 11.14, 9.33]
r = [.75, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]

scatter(r, velocity)
scatter(r, width, marker = '+')
print polyfit(log(r), log(velocity), 1)
print polyfit(log(r), log(width), 1)
show()

velocity2 = [.213, 8.436, 12.85, 17.5, 21.6, 25.5, 33.3, 40.45]
width2 = [4.9, 56, 57, 57, 58, 57.55, 57.824, 61.65]
velocity3 = [3.181, 5.572, 7.701, 15.4, 13.83, 12.2, 11.42]
width3 = [27.24, 28.66, 29.03, 34, 31.68, 29.75, 27.2]
c2 = [.01, .2, .3, .4, .5, .6, .8, 1.0] #r = .5
c3 = [.30, .50, .70, .85, .90, .95, 1.0] #r = 1.0. I must have made a mistake here. values for .85 and above are actually r values, not c values.

scatter(c2, velocity2)
scatter(c3, velocity3, marker='+')
show()
