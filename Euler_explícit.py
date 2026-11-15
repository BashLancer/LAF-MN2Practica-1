#importem llibreries
import numpy as np
import matplotlib.pyplot as plt

#definim els paràmetres
N=101
x=np.linspace(0,1,N) #creem N punts equiespaiats de 0 a 1
dx=1/(N-1) #Calculem les dimensions dels intervals de x
ta=0.025 #l'instant de temps final normalitzat.
gamma=0.51  #és la gamma=dt^/(dx^)^2
#Aïllant
dt=gamma*dx**2 #és llavors el temps per passa és aquest (és així per definició, per com ho hem discretitzat)
passos=round(ta/dt) #Dividim el temps total entre el temps que triga en donar cada passa per veure el num. de passos.
dt=ta/passos #amb aquest dt l'última passa no cau a ta (perque els passos no eren sencers) per tant, calculem la gamma que necessitariem pequè quadrés tot
gamma =dt/dx**2 #la nova gamma
