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


U_EE = Problema_Cauchy(Euler_Explicito, Oscilador, U0, t)
U_EI = Problema_Cauchy(Euler_Implicito, Oscilador, U0, t)
U_LP = Problema_Cauchy_LP(Leap_Frog, Oscilador, U0, t)
U_CN = Problema_Cauchy(Crank_Nicolson, Oscilador, U0, t)
U_RK4 = Problema_Cauchy(RK4, Oscilador, U0, t)




