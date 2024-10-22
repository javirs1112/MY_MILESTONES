from numpy import concatenate
from numpy.linalg import norm
from scipy.optimize import newton


def Problema_Cauchy(Esquema, F, U, N, t):
    
    for n in range(N):

        dt = t[n+1] - t[n]
        
        U[n+1,:] = Esquema(U[n,:], F, dt, n)

    return U

def Kepler(U,t):

    F = concatenate((U[2:4],-U[0:2]/norm(U[0:2])**3), axis = 0)

    return F
def Euler_Explicito(U, F, dt, t):

    return U + dt*F(U,t)

def RK4(U, F, dt, t):

    k1 = F(U,t)
    k2 = F(U + dt*k1/2, t + dt/2)
    k3 = F(U + dt*k2/2, t + dt/2)
    k4 = F(U + dt*k3, t + dt)

    return U + dt*(k1 + 2*k2 + 2*k3 + k4)/6

def Euler_Implicito(U, F, dt, t):

    def G(X):

        return X - U - dt*F(X, t)
    
    return newton(G, U) #Nota: Escribir Newton hecho por mí.

def Crank_Nicolson(U, F, dt, t):

    def G(X):

        return X - U - dt*(F(X, t)+F(U,t))
    
    return newton(G, U)