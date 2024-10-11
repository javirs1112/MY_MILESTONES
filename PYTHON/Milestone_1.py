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

t = linspace(t0,tf,N+1)

U = zeros([N+1,4])
U[0,:] = [x0, y0, vx0, vy0]

esquema = input("Introduce el esquema numérico deseado(Euler o RK4): ")

if esquema == "Euler": #ESQUEMA EULER EXPLÍCITO

    for n in range(N):

        U[n+1, 0] = U[n, 0] + deltat*U[n,2]
        U[n+1, 1] = U[n, 1] + deltat*U[n,3]
        U[n+1, 2] = U[n, 2] + deltat*(-U[n,0]/(U[n,0]**2 + U[n,1]**2)**(3/2))
        U[n+1, 3] = U[n, 3] + deltat*(-U[n,1]/(U[n,0]**2 + U[n,1]**2)**(3/2))

elif esquema == "RK4": #ESQUEMA RUNGE-KUTTA ORDEN 4

    for n in range(N):

       #Cálculo de k1 para cada variable de estado
        k1_x = U[n,2]
        k1_y = U[n,3]
        k1_Vx = (-U[n,0]/(U[n,0]**2 + U[n,1]**2)**(3/2))
        k1_Vy = (-U[n,1]/(U[n,0]**2 + U[n,1]**2)**(3/2))
       #Cálculo de k2 para cada variable de estado 
        k2_x = (U[n,2] + 0.5*deltat*k1_Vx)
        k2_y = (U[n,3] + 0.5*deltat*k1_Vy)
        k2_Vx = (-(U[n,0] + 0.5*deltat*k1_x)/((U[n,0] + 0.5*deltat*k1_x)**2 + (U[n,1] + 0.5*deltat*k1_y)**2)**(3/2))
        k2_Vy = (-(U[n,1] + 0.5*deltat*k1_y)/((U[n,0] + 0.5*deltat*k1_x)**2 + (U[n,1] + 0.5*deltat*k1_y)**2)**(3/2))
       #Cálculo de k3 para cada variable de estado 
        k3_x = (U[n,2] + 0.5*deltat*k2_Vx)
        k3_y = (U[n,3] + 0.5*deltat*k2_Vy)
        k3_Vx = (-(U[n,0] + 0.5*deltat*k2_x)/((U[n,0] + 0.5*deltat*k2_x)**2 + (U[n,1] + 0.5*deltat*k2_y)**2)**(3/2))
        k3_Vy = (-(U[n,1] + 0.5*deltat*k2_y)/((U[n,0] + 0.5*deltat*k2_x)**2 + (U[n,1] + 0.5*deltat*k2_y)**2)**(3/2))
       #Cálculo de k4 para cada variable de estado 
        k4_x = (U[n,2] + deltat*k3_Vx)
        k4_y = (U[n,3] + deltat*k3_Vy)
        k4_Vx = (-(U[n,0] + deltat*k3_x)/((U[n,0] + deltat*k3_x)**2 + (U[n,1] + deltat*k3_y)**2)**(3/2))
        k4_Vy = (-(U[n,1] + deltat*k3_y)/((U[n,0] + deltat*k3_x)**2 + (U[n,1] + deltat*k3_y)**2)**(3/2))

    #Calculamos el siguiente paso temporal de las variables de estado 

        U[n+1, 0] = U[n, 0] + (deltat/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x)
        U[n+1, 1] = U[n, 1] + (deltat/6)*(k1_y + 2*k2_y + 2*k3_y + k4_y)
        U[n+1, 2] = U[n, 2] + (deltat/6)*(k1_Vx + 2*k2_Vx + 2*k3_Vx + k4_Vx)
        U[n+1, 3] = U[n, 3] + (deltat/6)*(k1_Vy + 2*k2_Vy + 2*k3_Vy + k4_Vy)

else:
    print("Esquema no válido")


plt.plot(U[:,0],U[:,1])
plt.show()


