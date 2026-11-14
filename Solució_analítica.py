# S'importen llibreries:
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

# S'afegeixen les dades del problema:
d = 0.02                              # Distància/gruix (m)
V = 40                                # Voltatge (v)
c_v = 3698                            # Calor específic (J/(Kg*K))
p = 1081                              # Densitat (kg/m^3)
k = 0.56                              # Conductivitat térmica
sigma = 0.472                         # Conductividad eléctrica
N = 101                               # Número de punts, mallat espaial

# Es calculen els paràmetres:
P_ext = (sigma*(V**2))/(2*d**2)       # Fonts externes de calor del sistema
t_0 = ((d**2)*c_v*p)/k                # Temps normalitzat
T_0 = ((d**2)*P_ext)/k                # Temperatura normalitzada

# Es normalitza les variables espaials i temporals:
d_norm = 0.02/d                       # Distància normalitzada
dx = d_norm/100                       # Variació del mallat espaial
ta = 0.025                            # Valor del mallat temporal final

# Es ressol analíticament el problema: 

# Es defineix la temperatura inicial y es crea un mallat espaial:
x = np.linspace(0, 2, 101)
T_inicial = 36.5 / T_0
h = []

# Es calcula la distribució de temperatures, mitjançant una série de Fourier:
for a in range(0, 101):
  sumatorio = 0
  for i in range(1,9999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  h.append(T_inicial + sumatorio)
  
# Es recupera la temperatura:
l = np.array(h) * T_0

# Es representa gràficament:
plt.plot(x,l)
plt.xlabel('Posición (cm)')
plt.ylabel('Temperatura (ºC)')
plt.title('Solución analítica')
plt.axvspan(0.75, 1.25, facecolor ='lightblue', alpha=0.3)
plt.axvspan(0, 0.75, facecolor ='lightcoral', alpha=0.3)
plt.axvspan(1.25, 2., facecolor ='lightcoral', alpha=0.3)
plt.show()
