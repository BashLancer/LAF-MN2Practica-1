# S'afegeixen les dades del problema:
d = 0.02                              # Distància/gruix (m)
V = 40                                # Voltatge (v)
C_v = 3698                            # Calor específic (J/(Kg*K))
p = 1081                              # Densitat (kg/m^3)
k = 0.56                              # Conductivitat térmica
sigma = 0.472                         # Conductividad eléctrica
N = 101                               # Número de punts, mallat espaial

P_ext = (sigma*(V**2))/(2*d**2)       # Fonts externes de calor del sistema
t_0 = ((d**2)*c_v*p)/k                # Temps normalitzat
T_0 = ((d**2)*P_ext)/k                # Temperatura normalitzada

d_norm = 0.02/d                       # Distància normalitzada
dx = d_norm/100                       # Variació del mallat espaial
ta = 0.025                            # Valor del mallat temporal final
