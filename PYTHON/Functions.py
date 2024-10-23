from numpy import concatenate, zeros, linspace
from numpy.linalg import norm
from scipy.optimize import newton


def Problema_Cauchy(Esquema, F, U0, t):
    
    N = len(t) - 1
    U = zeros([N+1, len(U0)])
    U[0,:] = U0

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

def Richardson(U0, F, Problema, Esquema, t, q):  #Este método esta escrito sólo para usarse con dt1 = dt2/2, generalizar en el futuro

    t1 = t
    t2 = linspace(t[0], t[-1], 2*len(t1))

    Error = zeros([len(t1), len(U0)])           

    U1 = Problema(Esquema, F, U0, t1)
    U2 = Problema(Esquema, F, U0, t2)

    for i in range(len(t)):
        Error[i,:] = (U2[2*i,:] - U1[i,:]) / (1 - 1/2**q)

    return Error

