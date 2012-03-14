from multiprocessing import Pool

def f(x):
    return x*x

def y : return x*100

pool = Pool(processes=4)              # start 4 worker processes
print pool.map(y, range(1000))

