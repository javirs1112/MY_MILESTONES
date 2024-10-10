from numpy import array, zeros, linspace
from numpy.linalg import norm
import matplotlib.pyplot as plt

#DATOS Y CONDICIONES INICIALES

x0 = 1
y0 = 0
vx0 = 0
vy0 = 1

t0 = 0
tf = 20
N = 200
deltat = (tf-t0)/N

t = linspace(t0,tf,N)


