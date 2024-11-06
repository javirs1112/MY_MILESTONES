from numpy import array, zeros, linspace
from numpy.linalg import norm
import matplotlib.pyplot as plt
from Functions import Problema_Cauchy, Kepler, Oscilador, Euler_Explicito, RK4, Euler_Implicito, Crank_Nicolson

#DATOS Y CONDICIONES INICIALES

U0 = array([1, 0])


t0 = 0
tf = 10
N = 1000

t = linspace(t0,tf,N+1)

U = zeros([N+1,len(U0)])
U[0,:] = U0

U = Problema_Cauchy(Euler_Implicito, Oscilador, U0, t)

plt.axis('equal') 
plt.plot(t, U[:,0])
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