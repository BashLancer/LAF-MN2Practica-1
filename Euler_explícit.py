#importem llibreries
import numpy as np
import matplotlib.pyplot as plt

#definim els paràmetres
N=101
x=np.linspace(0,1,N) #creem N punts equiespaiats de 0 a 1
delta_x=1/(N-1) #Calculem les dimensions dels intervals de x
t_a=0.025 #l'instant de temps final normalitzat.
