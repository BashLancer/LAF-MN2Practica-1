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


















# SOLUCIÓ NUMÈRICA: EULER IMPLICIT

# Euler implícit amb bucle de Gauss-Seidel per dx ** 2:

# S'escriu les condicions inicials:
dt = dx ** 2
T_inicial = 36.5 / T_0
iteraciones = int(ta / dt) + 1

# Es crea una matriu de temperatures i s'aplica les condicions inicials i de contorn:
Matriz_4 = np.ones((iteraciones, N))
for i in range(N):
  Matriz_4[0, i] = T_inicial
for j in range(iteraciones):
  Matriz_4[j, 0] = T_inicial
  Matriz_4[j, -1] = T_inicial
  
M = Matriz_4.copy()

# Es ressol numèricament mitjançant Gauss-Seidel:
for k in range(1000):
  for fila in range(1, iteraciones):
    for columna in range(1, N-1):
      M[fila, columna] = 1 / (1 + 2 * dt /(dx ** 2)) * (((M[fila, columna + 1]) + (M[fila, columna - 1])) * (dt / (dx ** 2)) + dt + (M[fila - 1, columna]))
    Matriz_4 = M.copy()
    
# Es calcula la solució analítica:
ta = (iteraciones - 1) * dt
T_inicial = 36.5 / T_0
h= []
for a in range(0, 101):
  sumatorio = 0
  for i in range(1,999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  h.append(T_inicial + sumatorio)
l_4 = np.array(h) * T_0

# Es calcula els errors i es compara gràficament:
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.set_xlabel('Posición (cm)')
ax.set_ylabel('Temperatura (ºC)')
plt.plot(x, l_4, label = 'Solución analítica')
plt.plot(x, Matriz_4[-1,:] * T_0, label = 'Gauss_Seidel dt = dx ** 2')
ax.set_title('Solución analítica vs Solución Gauss_Seidel')
ax.axvspan(0.75, 1.25, facecolor ='lightblue', alpha=0.5, label='Región 1')
ax.axvspan(0, 0.75, facecolor ='lightcoral', alpha=0.5, label='Región 2')
ax.axvspan(1.25, 2., facecolor ='lightcoral', alpha=0.5)
ax.legend()
plt.show()

error_4 = np.abs(Matriz_4[-1,:] * T_0 - l_4)
plt.plot(x, error_4)
plt.title('Error en función de la posición')
plt.xlabel('Posición (cm)')
plt.ylabel('Error (ºC)')
plt.show()

# Euler implícit amb bucle de Gauss-Seidel per 0.5 * dx ** 2:

# S'escriu les condicions inicials:
dt = 0.5 * dx ** 2
T_inicial = 36.5 / T_0
iteraciones = int(ta / dt) + 1

# Es crea una matriu de temperatures i s'aplica les condicions inicials i de contorn:
Matriz_5 = np.ones((iteraciones, N))
for i in range(N):
  Matriz_5[0, i] = T_inicial
for j in range(iteraciones):
  Matriz_5[j, 0] = T_inicial
  Matriz_5[j, -1] = T_inicial

M = Matriz_5.copy()

# Es ressol numèricament mitjançant Gauss-Seidel:
for k in range(1000):
  for fila in range(1, iteraciones):
    for columna in range(1, N-1):
        M[fila, columna] = 1 / (1 + 2 * dt /(dx ** 2)) * (((M[fila, columna + 1]) + (M[fila, columna - 1])) * (dt / (dx ** 2)) + dt + (M[fila - 1, columna]))
    Matriz_5 = M.copy()

# Es calcula la solució analítica:
x = np.linspace(0, 2, 101)
ta = (iteraciones - 1) * dt
T_inicial = 36.5 / T_0
h= []
for a in range(0, 101):
  sumatorio = 0
  for i in range(1,999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  h.append(T_inicial + sumatorio)
l_5 = np.array(h) * T_0

# Es calcula els errors i es compara gràficament:
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.set_xlabel('Posición (cm)')
ax.set_ylabel('Temperatura (ºC)')
plt.plot(x, l_5, label = 'Solución analítica')
plt.plot(x, Matriz_5[-1,:] * T_0, label = 'Gauss_Seidel dt = 0.5 dx ** 2')
ax.set_title('Solución analítica vs Solución Gauss_Seidel')
ax.axvspan(0.75, 1.25, facecolor ='lightblue', alpha=0.5, label='Región 1')
ax.axvspan(0, 0.75, facecolor ='lightcoral', alpha=0.5, label='Región 2')
ax.axvspan(1.25, 2., facecolor ='lightcoral', alpha=0.5)
ax.legend()
plt.show()

error_5 = np.abs(Matriz_5[-1,:] * T_0 - l_5)
plt.plot(x, error_5)
plt.title('Error en función de la posición')
plt.xlabel('Posición (cm)')
plt.ylabel('Error (ºC)')
plt.show()
