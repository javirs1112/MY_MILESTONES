from numpy import array, zeros, linspace
from numpy.linalg import norm
import matplotlib.pyplot as plt
from Functions import Problema_Cauchy, Kepler, Oscilador, Euler_Explicito, RK4, Euler_Implicito, Crank_Nicolson, Leap_Frog, Problema_Cauchy_LP

#DATOS Y CONDICIONES INICIALES

U0 = array([1, 0])


t0 = 0
tf = 10
N = 1000

t = linspace(t0,tf,N+1)

U = zeros([N+1,len(U0)])
U[0,:] = U0


U = Problema_Cauchy(Leap_Frog, Oscilador, U0, t)

plt.plot(t, U[:,0])
plt.show()





