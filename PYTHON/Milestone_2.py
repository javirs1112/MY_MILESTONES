from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt

def Problema_Cauchy(Esquema, F, U, N, dt, t):

    for n in range(N):
        
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

#DATOS Y CONDICIONES INICIALES

U0 = array([1, 0, 0, 1])


t0 = 0
tf = 20
N = 200
dt = (tf-t0)/N

t = linspace(t0,tf,N+1)

U = zeros([N+1,4])
U[0,:] = U0

U = Problema_Cauchy(Euler_Explicito, Kepler, U, N, dt, t)

plt.plot(U[:,0],U[:,1])
plt.show()

# for n in range(N):

#     U[n+1, 0] = U[n, 0] + dt*U[n,2]
#     U[n+1, 1] = U[n, 1] + dt*U[n,3]
#     U[n+1, 2] = U[n, 2] + dt*(-U[n,0]/(U[n,0]**2 + U[n,1]**2)**(3/2))
#     U[n+1, 3] = U[n, 3] + dt*(-U[n,1]/(U[n,0]**2 + U[n,1]**2)**(3/2))



# for n in range(N):

#     #Cálculo de k1 para cada variable de estado
#     k1_x = U[n,2]
#     k1_y = U[n,3]
#     k1_Vx = (-U[n,0]/(U[n,0]**2 + U[n,1]**2)**(3/2))
#     k1_Vy = (-U[n,1]/(U[n,0]**2 + U[n,1]**2)**(3/2))
#     #Cálculo de k2 para cada variable de estado 
#     k2_x = (U[n,2] + 0.5*dt*k1_Vx)
#     k2_y = (U[n,3] + 0.5*dt*k1_Vy)
#     k2_Vx = (-(U[n,0] + 0.5*dt*k1_x)/((U[n,0] + 0.5*dt*k1_x)**2 + (U[n,1] + 0.5*dt*k1_y)**2)**(3/2))
#     k2_Vy = (-(U[n,1] + 0.5*dt*k1_y)/((U[n,0] + 0.5*dt*k1_x)**2 + (U[n,1] + 0.5*dt*k1_y)**2)**(3/2))
#     #Cálculo de k3 para cada variable de estado 
#     k3_x = (U[n,2] + 0.5*dt*k2_Vx)
#     k3_y = (U[n,3] + 0.5*dt*k2_Vy)
#     k3_Vx = (-(U[n,0] + 0.5*dt*k2_x)/((U[n,0] + 0.5*dt*k2_x)**2 + (U[n,1] + 0.5*dt*k2_y)**2)**(3/2))
#     k3_Vy = (-(U[n,1] + 0.5*dt*k2_y)/((U[n,0] + 0.5*dt*k2_x)**2 + (U[n,1] + 0.5*dt*k2_y)**2)**(3/2))
#     #Cálculo de k4 para cada variable de estado 
#     k4_x = (U[n,2] + dt*k3_Vx)
#     k4_y = (U[n,3] + dt*k3_Vy)
#     k4_Vx = (-(U[n,0] + dt*k3_x)/((U[n,0] + dt*k3_x)**2 + (U[n,1] + dt*k3_y)**2)**(3/2))
#     k4_Vy = (-(U[n,1] + dt*k3_y)/((U[n,0] + dt*k3_x)**2 + (U[n,1] + dt*k3_y)**2)**(3/2))

# #Calculamos el siguiente paso temporal de las variables de estado 

#     U[n+1, ] = U[n, 0] + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x)