#importem llibreries
import numpy as np
import matplotlib.pyplot as plt

N=101
x=np.linspace(0,1,N) #creem N punts equiespaiats de 0 a 1
dx=1/(N-1) #Calculem les dimensions dels intervals de x
ta=0.025 #l'instant de temps final normalitzat.

def solucio(gamma):


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
      T[j]=T_anterior[j]+gamma*(T_anterior[j+1]-2*T_anterior[j]+T_anterior[j-1])+dt
  return T
#ara trobem l'error numèric absolut (no s'especifica quin voleu)
#Per fer-ho necessitem la solució analítica

T_anal=np.zeros(N)

for a in range(0, 101):
  sumatorio = 0
  for i in range(1,9999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  T_anal[a] = sumatorio #el meu company havia creat una llista però jo necessito arrays.

#error de la gamma 1
error_abs1= np.abs(T_anal-solucio(0.51))

#error de la gamma 2
error_abs2= np.abs(T_anal-solucio(0.49))

#error de la gamma 3
error_abs3= np.abs(T_anal-solucio(0.25))


#grafiquem. Separem gamma=0,51 de les altres dues perque sinó, no es veu bé.

#Fem la part esquerra de la figura
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,4))
ax1.plot(x,T_anal,color='black',linestyle='-',label='S.Analítica')
ax1.plot(x,solucio(0.49),color='red',linestyle='--',label=r'$\gamma = 0.49$')
ax1.plot(x,solucio(0.25),color='blue',linestyle=':',label=r'$\gamma = 0.25$')
ax1.set_xlabel(r'$\hat{x}$')
ax1.set_ylabel(r'$\hat{T}$')
ax1.set_title(r'Solució a $\hat{t}_a = 0.025$')
ax1.legend()
#Fem ara la part dreta de la figura, l'error.
ax2.plot(x,error_abs2,color='red',label=r'$\gamma = 0.49$')
ax2.plot(x,error_abs3,color='blue',label=r'$\gamma = 0.25$')
ax2.set_xlabel(r'$\hat{x}$')
ax2.set_ylabel('Error absolut')
ax2.set_title('Error numèric')
ax2.legend()



#part esquerra de gamma=0,51

fig,(ax1,ax2) =plt.subplots(1,2,figsize=(10,4))
ax1.plot(x,T_anal,color='black',label='S.Analítica')
ax1.plot(x,T_anal-solucio(0.51),color='green',label=r'$\gamma = 0.51$')
ax1.set_xlabel(r'$\hat{x}$')
ax1.set_ylabel(r'$\hat{T}$')
ax1.set_title(r'Solució a $\hat{t}_a = 0.025$')
ax1.legend()
#part dreta gamma=0,51
ax2.plot(x, error_abs1,color='green')
ax2.set_xlabel(r'$\hat{x}$')
ax2.set_ylabel('Error absolut')
ax2.set_title('Error numèric')
plt.show()

