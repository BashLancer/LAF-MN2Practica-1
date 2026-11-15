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

# Es compara la representació gràfica de la solució numèrica amb els mallats corresponents, en funció de la solució analítica:
plt.plot(x, Matriz_5[-1,:] * T_0, label = 'dt = 0.5 dx**2')
plt.plot(x, Matriz_4[-1,:] * T_0, label = 'dt = dx**2')
plt.plot(x, l_5, label = 'Sol analítica dt = 0.5 dx**2')
plt.plot(x, l_4, label = 'Sol analítica dt = dx**2')
plt.axvspan(0, 0.75, alpha=0.5, color='lightcoral')
plt.axvspan(0.75, 1.25, alpha=0.5, color='lightblue')
plt.axvspan(1.25, 2, alpha=0.5, color='lightcoral')
plt.xlabel('Posición (cm)')
plt.ylabel('Temperatura (ºC)')
plt.legend()
plt.show()

plt.plot(x,error_4, label = 'Error 1')
plt.plot(x, error_5, label = 'Error 0.5')
plt.xlabel('Posición (cm)')
plt.ylabel('Error')
plt.legend()
plt.show()

# COMPARACIÓ AMB EL MÉTODE D'EULER EXPLÍCIT:

# Euler implícit amb bucle de Gauss-Seidel per 0.51 * dx ** 2:

# S'escriu les condicions inicials:
dt = 0.51 * dx ** 2
T_inicial = 36.5 / T_0
iteraciones = int(ta / dt) + 1

# Es crea una matriu de temperatures i s'aplica les condicions inicials i de contorn:
Matriz_6 = np.ones((iteraciones, N))

for i in range(N):
  Matriz_6[0, i] = T_inicial
for j in range(iteraciones):
  Matriz_6[j, 0] = T_inicial
  Matriz_6[j, -1] = T_inicial
  
M = Matriz_6.copy()

# Es ressol numèricament mitjançant Gauss-Seidel:
for k in range(1000):
  for fila in range(1, iteraciones):
    for columna in range(1, N-1):
      M[fila, columna] = 1 / (1 + 2 * dt /(dx ** 2)) * (((M[fila, columna + 1]) + (M[fila, columna - 1])) * (dt / (dx ** 2)) + dt + (M[fila - 1, columna]))
    Matriz_6 = M.copy()

# Es calcula la solució analítica:
ta = (iteraciones - 1) * dt
T_inicial = 36.5 / T_0
h= []
for a in range(0, 101):
  sumatorio = 0
  for i in range(1,999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  h.append(T_inicial + sumatorio)
l_6 = np.array(h) * T_0

# Es calcula els errors i es compara gràficament:
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.set_xlabel('Posición (cm)')
ax.set_ylabel('Temperatura (ºC)')
plt.plot(x, l_6, label = 'Solución analítica')
plt.plot(x, Matriz_6[-1,:] * T_0, label = 'Gauss_Seidel dt = 0.51 dx ** 2')
ax.set_title('Solución analítica vs Solución Gauss_Seidel')
ax.axvspan(0.75, 1.25, facecolor ='lightblue', alpha=0.5, label='Región 1')
ax.axvspan(0, 0.75, facecolor ='lightcoral', alpha=0.5, label='Región 2')
ax.axvspan(1.25, 2., facecolor ='lightcoral', alpha=0.5)
ax.legend()
plt.show()

error_6 = np.abs(Matriz_6[-1,:] * T_0 - l_6)
plt.plot(x, error_6)
plt.title('Error en función de la posición')
plt.xlabel('Posición (cm)')
plt.ylabel('Error (ºC)')
plt.show()

# Euler implícit amb bucle de Gauss-Seidel per 0.49 * dx ** 2:

# S'escriu les condicions inicials:
dt = 0.49 * dx ** 2
T_inicial = 36.5 / T_0
iteraciones = int(ta / dt) + 1

# Es crea una matriu de temperatures i s'aplica les condicions inicials i de contorn:
Matriz_7 = np.ones((iteraciones, N))

for i in range(N):
  Matriz_7[0, i] = T_inicial
for j in range(iteraciones):
  Matriz_7[j, 0] = T_inicial
  Matriz_7[j, -1] = T_inicial

M = Matriz_7.copy()

# Es ressol numèricament mitjançant Gauss-Seidel:
for k in range(1000):
  for fila in range(1, iteraciones):
    for columna in range(1, N-1):
        M[fila, columna] = 1 / (1 + 2 * dt /(dx ** 2)) * (((M[fila, columna + 1]) + (M[fila, columna - 1])) * (dt / (dx ** 2)) + dt + (M[fila - 1, columna]))
    Matriz_7 = M.copy()
    
# Es calcula la solució analítica:
ta = (iteraciones - 1) * dt
T_inicial = 36.5 / T_0
h= []
for a in range(0, 101):
  sumatorio = 0
  for i in range(1,999):
    sumatorio += (4 * (1 - np.exp(- np.pi ** 2 * ((2 * i - 1) ** 2) * ta))/(((2 * i - 1) ** 3) * np.pi ** 3)) * np.sin((2 * i - 1) * np.pi * a * dx)
  h.append(T_inicial + sumatorio)
l_7 = np.array(h) * T_0

# Es calcula els errors i es compara gràficament:
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.set_xlabel('Posición (cm)')
ax.set_ylabel('Temperatura (ºC)')
plt.plot(x, l_7, label = 'Solución analítica')
plt.plot(x, Matriz_7[-1,:] * T_0, label = 'Gauss_Seidel dt = 0.49 dx ** 2')
ax.set_title('Solución analítica vs Solución Gauss_Seidel')
ax.axvspan(0.75, 1.25, facecolor ='lightblue', alpha=0.5, label='Región 1')
ax.axvspan(0, 0.75, facecolor ='lightcoral', alpha=0.5, label='Región 2')
ax.axvspan(1.25, 2., facecolor ='lightcoral', alpha=0.5)
ax.legend()
plt.show()

error_7 = np.abs(Matriz_7[-1,:] * T_0 - l_7)
plt.plot(x, error_7)
plt.title('Error en función de la posición')
plt.xlabel('Posición (cm)')
plt.ylabel('Error (ºC)')
plt.show()
