from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt
from Milestone_2 import Kepler, Euler_Explicito, Problema_Cauchy
from scipy.optimize import newton

def Euler_Implicito(U, F, dt, t):

    def G(X):

        return X - U - dt*F(X, t)
    
    return newton(G, U) #Nota: Escribir Newton hecho por mí.

def Crank_Nicolson(U, F, dt, t):

    def G(X):

        return X - U - dt*(F(X, t)+F(U,t))
    
    return newton(G, U)

#DATOS Y CONDICIONES INICIALES

U0 = array([1, 0, 0, 1])


t0 = 0
tf = 7
N = 800
dt = (tf-t0)/N

t = linspace(t0,tf,N+1)

U = zeros([N+1,4])
U[0,:] = U0

U = Problema_Cauchy(Euler_Implicito, Kepler, U, N, t)


plt.plot(U[:,0], U[:,1])
plt.show()
