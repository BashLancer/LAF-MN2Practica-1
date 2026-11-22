#HE REUTILITZAT PART DEL CODI D'EULER EXPLÍCIT

#importem llibreries
import numpy as np
import matplotlib.pyplot as plt

#Afegim els paràmetres del sistema
L=0.02
V=40
cv=3698
rho=1081
k=0.56
sigma=0.472
Tc=36.5
#Calulem P_ext
P_ext = sigma*V**2/(2*L**2)
#i també els temps i temperatures que normalitzen l'equació de calor
t_0 = cv*rho*L**2/k
T_0 = P_ext*L**2/k
#afegim ara les temperatures que diu el guió que no pot sobrepasar cada tipus de teixit
T50=(50-Tc)/T_0
T80=(80-Tc)/T_0

N=101
x=np.linspace(0,1,N) #creem N punts equiespaiats de 0 a 1
dx=1/(N-1) #Calculem les dimensions dels intervals de x
ta=0.1 #Ara necessitem un temps gran per veure quan es compleixen les condicions.

gamma=0.25


#Aïllant
dt=gamma*dx**2 #és llavors el temps per passa és aquest (és així per definició, per com ho hem discretitzat)
passos=round(ta/dt) #Dividim el temps total entre el temps que triga en donar cada passa per veure el num. de passos.
dt=ta/passos #amb gamma 0.51 l'última passa no cau a ta, per tant, calculer la gamma que necessitariem pequè quadrés tot
gamma =dt/dx**2 #gamma real

#Tenim 101 temperatures i les hem de guardar en algun lloc, aixó que creem un array
T=np.zeros(N)
Totes_temperatures=np.zeros((passos+1,N)) #necessitem totes, no com abans que necessitavem solament t_a. Passos+1 files i N columnes

  #Utilitzem un bucle per calcular la nova temperatura a cada punt un instant després amb la formulilla

for i in range(passos):
  T_anterior = T.copy()
  for j in range(1,len(T)-1):
    T[j]=T_anterior[j]+gamma*(T_anterior[j+1]-2*T_anterior[j]+T_anterior[j-1])+dt
  Totes_temperatures[i]=T
    #creem un array on ficar tots els temps(en segons) fins a ta i desfent la normalització
temps=np.zeros(passos+1)
for i in range(passos+1):
  temps[i]=i*dt*t_0


#Però ara hem de saber quin és l'instant de temps on es compleixen les condicions de temperatura per matar el teixit defectuós i no fer mal al sa.
T_malalt=(x>=0.375)&(x<=0.625) #posició de cada tipus de teixit
T_sa=(x<0.375)|(x>0.625) #O una cosa o l'altre

frac_malalt=np.mean(Totes_temperatures[:,T_malalt]>=T50,axis=1) #Quina part del teixit malalt ha arribat ja als 50 graus a cada instant?
frac_sa=np.mean(Totes_temperatures[:,T_sa]>=T50,axis=1) #axis=1 fa la mesura de la mitjana al llarg de cada fila

i_malalt=np.argmax(frac_malalt>= 1) #Retorna la posició del primer true on TOT el teixit malalt ha arribat a 50 graus
i_sa=np.argmax(frac_sa>0) #Retorna la posició del primer true on EL PRIMER punt de teixit sa ha arribat a 50 graus

print('Tot el malalt >50 ºC a', temps[i_malalt], 's')
print('Primer punt sa > 50 ºC a', temps[i_sa], 's')

#grafiquem
plt.figure(figsize=(10, 6))
plt.plot(x*2,Tc+T_0*Totes_temperatures[i_sa])
plt.xlabel('Posició (cm)')
plt.ylabel('Temperatura (ºC)')
plt.show()
