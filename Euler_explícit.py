#importem llibreries
import numpy as np
import matplotlib.pyplot as plt

N=101
x=np.linspace(0,1,N) #creem N punts equiespaiats de 0 a 1
dx=1/(N-1) #Calculem les dimensions dels intervals de x
ta=0.025 #l'instant de temps final normalitzat.

gamma=0.51  #és la gamma=dt^/(dx^)^2
#Aïllant
dt=gamma*dx**2 #és llavors el temps per passa és aquest (és així per definició, per com ho hem discretitzat)
passos=round(ta/dt) #Dividim el temps total entre el temps que triga en donar cada passa per veure el num. de passos.
dt=ta/passos #amb gamma 0.51 l'última passa no cau a ta, per tant, calculer la gamma que necessitariem pequè quadrés tot
gamma =dt/dx**2 #gamma real

#Tenim 101 temperatures i les hem de guardar en algun lloc, aixó que creem un array
T=np.zeros(N)

#Utilitzem un bucle per calcular la nova temperatura a cada punt un instant després amb la formulilla

for i in range(passos):
  T_anterior = T.copy()
  for j in range(1,len(T)-1):
    T[j]=T_anterior[j]+gamma*(T[j+1]-2*T[j]+T[j-1])+dt

#ara trobem l'error numèric absolut (no s'especifica quin voleu)
#Per fer-ho necessitem la solució analítica

T_anal=np.zeros(N)

for a in range(0, 101):
  sumatorio = 0
  for i in range(1,9999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  T_anal[a] = sumatorio #el meu company havia creat una llista però jo necessito arrays

error_abs= np.abs(T_anal-T)

