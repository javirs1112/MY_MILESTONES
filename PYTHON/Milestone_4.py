from numpy import array, zeros, linspace, shape
from numpy.linalg import norm
import matplotlib.pyplot as plt
from Functions import Problema_Cauchy, Kepler, Oscilador, Euler_Explicito, RK4, Euler_Implicito, Crank_Nicolson, Leap_Frog, Region_Estabilidad

#DATOS Y CONDICIONES INICIALES

U0 = array([1, 0])


t0 = 0
tf = 20
N = 100

t = linspace(t0,tf,N+1)

U = zeros([N+1,len(U0)])
U[0,:] = U0


# U = Problema_Cauchy(Leap_Frog, Oscilador, U0, t)

# plt.plot(t, U[:,0])
# plt.show()

x, y, rho = Region_Estabilidad(Leap_Frog, N, -4, 2, -4, 4)

print(shape(rho))
plt.contour(x, y, rho, linspace(0, 1, 11))
plt.grid()
plt.title('Región de estabilidad')
plt.xlabel('Im(w)')
plt.ylabel('Re(w)')
plt.show()





