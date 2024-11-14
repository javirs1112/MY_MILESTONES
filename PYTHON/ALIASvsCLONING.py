from numpy import array
V = array([1,2,3])
pV = V ## ALIAS: sirve para dar dos nombres al mismo objeto

pV[0] = 4

print(V)
print(id(pV)) # Mismo id: dos objetos o pointers que apuntan al mismo espacio de memoria (mismo vector)
print(id(V))

U = V.copy() ## COPY: sirve para crear un duplicado de V
print(id(U)) # Diferente id que V