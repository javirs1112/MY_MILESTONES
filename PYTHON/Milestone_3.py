from Functions import Problema_Cauchy, Kepler, Euler_Explicito, RK4, Euler_Implicito, Crank_Nicolson, Richardson
from numpy import array, linspace
import matplotlib.pyplot as plt

U0 = array([1, 0, 0, 1])
t0 = 0
tf = 10
N = 1000
t = linspace(t0, tf, 1000)

input = input("Introduce el esquema numérico deseado(EE, RK4, EI, CN): ")

if input == "EE":
    esquema = Euler_Explicito
    q = 1

elif input == "RK4":
    esquema = RK4
    q = 4

elif input == "EI":
    esquema = Euler_Implicito
    q = 4

elif input == "CN":
    esquema = Crank_Nicolson
    q = 1
else:
    print("Esquema no válido")

Error = Richardson(U0, Kepler, Problema_Cauchy, esquema, t, q)

plt.axis('equal') 
plt.xlabel('tiempo')
plt.ylabel('Error')
plt.plot(t, Error[:,0], '-b', label = 'X')
plt.plot(t, Error[:,1], '-r', label = 'Y')
plt.plot(t, Error[:,2], '--m', label = 'Vx')
plt.plot(t, Error[:,3], '--c', label = 'Vy')
plt.legend()
plt.show()