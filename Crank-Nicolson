import numpy as np
import matplotlib.pyplot as plt



# Primer de tot definim els parámetres del problema

N = 101                  # Nombre de punts espacials
t_final = 0.025          # Temps final adimensional
L = 1.0                  # Longitud del domini normalitzat

dx = L / (N - 1)
x = np.linspace(0, L, N)

# Tolerancia i nombre màxim d'iteracions de Gauss-Seidel
tolerancia = 1e-12
max_iter = 100000

# Passos temporals que demana el guió
factors_dt = [1.0, 0.5]

# Ara calculem la solució analítica

def solucion_analitica(x, t, termes=399):
    """
    Calcula la solucio analitica de l'equacio:

        dT/dt = d2T/dx2 + 1

    amb condicions inicials i de contorn nul·les.
    """

    T = np.zeros_like(x, dtype=float)

    # Serie de Fourier: nomes apareixen els termes senars
    for n in range(1, termes + 1, 2):

        coeficient = 4.0 / (n * np.pi)

        T += (
            coeficient
            * (1.0 - np.exp(-n**2 * np.pi**2 * t))
            / (n**2 * np.pi**2)
            * np.sin(n * np.pi * x)
        )

    return T


# Començem amb el métode de Gaus-Seidel

def gauss_seidel(b, gamma, x_inicial):
    """
    Resol el sistema tridiagonal:

        -gamma/2 * x[j-1]
        + (1+gamma) * x[j]
        - gamma/2 * x[j+1] = b[j]

    # Mitjancant el metode iteratiu de Gauss-Seidel.

    Les condicions de contorn son x(0) = x(1) = 0.
    """

    x_sol = x_inicial.copy()

    # Nombre de punts interiors
    m = len(b)

    for iteracio in range(max_iter):

        x_anterior = x_sol.copy()

        # Recorrem els punts interiors
        for j in range(m):

            esquerra = 0.0 if j == 0 else x_sol[j - 1]
            dreta = 0.0 if j == m - 1 else x_sol[j + 1]

            x_sol[j] = (
                b[j] + (gamma / 2.0) * (esquerra + dreta)
            ) / (1.0 + gamma)

        # El criteri de convergencia serà:
        error_iteratiu = np.max(np.abs(x_sol - x_anterior))

        if error_iteratiu < tolerancia:
            return x_sol, iteracio + 1

    raise RuntimeError(
        "Gauss-Seidel no ha convergit. "
        "Cal augmentar max_iter o revisar la tolerancia."
    )

# Fem el métodode de Crank-Nicolson

def crank_nicolson(dt):

    gamma = dt / dx**2

    # Nombre de passos temporals
    Nt = round(t_final / dt)

    # Comprovem que arribem exactament al temps final
    assert np.isclose(Nt * dt, t_final)

    # Condicio inicial: temperatura nul·la normalitzada
    T = np.zeros(N)

    # Evolucio temporal
    temps = np.zeros(Nt + 1)

    # Guardem la temperatura final
    for n in range(Nt):

        # Terme independent del sistema

        b = (
            (gamma / 2.0) * T[:-2]
            + (1.0 - gamma) * T[1:-1]
            + (gamma / 2.0) * T[2:]
            + dt
        )

        # RESOLUCIÓ DEL SISTEMA AMB GAUSS-SEIDEL: 

        # Utilitzem la solucio de l'instant anterior com a aproximacio inicial
      
        T_inicial = T[1:-1].copy()

        T_interior, iteracions = gauss_seidel(
            b, gamma, T_inicial
        )

        # Actualitzem la distribucio de temperatures
        T_nova = np.zeros(N)
        T_nova[1:-1] = T_interior

        T = T_nova

        temps[n + 1] = (n + 1) * dt

    return T, temps, gamma


# Executem els dos casos

resultats = []

for factor in factors_dt:

    dt = factor * dx**2

    print("\n====================================")
    print("CRANK-NICOLSON")
    print("====================================")

    print(f"Factor temporal: {factor}")
    print(f"dx = {dx:.8f}")
    print(f"dt = {dt:.8f}")

    T_num, temps, gamma = crank_nicolson(dt)

    # Solucio analitica al temps final
    T_exacta = solucion_analitica(x, t_final)

    # Error absolut
    error = np.abs(T_num - T_exacta)

    error_maxim = np.max(error)
    error_mitja = np.mean(error)

    print(f"Gamma = {gamma:.4f}")
    print(f"Nombre de passos temporals = {len(temps)-1}")
    print(f"Temperatura maxima numerica = {np.max(T_num):.8f}")
    print(f"Error maxim = {error_maxim:.8e}")
    print(f"Error mitja = {error_mitja:.8e}")

    resultats.append({
        "factor": factor,
        "dt": dt,
        "T": T_num.copy(),
        "error": error.copy(),
        "gamma": gamma
    })

    # Grafiquem les temperatures
    plt.figure(figsize=(9, 5))

    plt.plot(
        x, T_num,
        label="Crank-Nicolson",
        linewidth=2
    )

    plt.plot(
        x, T_exacta,
        "--",
        label="Solucio analitica",
        linewidth=2
    )

    plt.xlabel(r"$\hat{x}$")
    plt.ylabel(r"$\hat{T}$")
    plt.title(f"Crank-Nicolson: dt = {dt:.6f}")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        f"CN_temperatura_{factor}.png",
        dpi=300
    )

    plt.show()

    # Grafiquem els errors 

    plt.figure(figsize=(9, 5))

    plt.plot(
        x, error,
        linewidth=2,
        label="Error absolut"
    )

    plt.xlabel(r"$\hat{x}$")
    plt.ylabel("Error absolut")
    plt.title(f"Error de Crank-Nicolson: dt = {dt:.6f}")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(
        f"CN_error_{factor}.png",
        dpi=300
    )

    plt.show()


# Comparem els dos casos temporals

plt.figure(figsize=(9, 5))

for resultat in resultats:

    plt.plot(
        x,
        resultat["T"],
        linewidth=2,
        label=f'dt = {resultat["dt"]:.6f}'
    )

plt.xlabel(r"$\hat{x}$")
plt.ylabel(r"$\hat{T}$")
plt.title("Comparacio dels dos casos de Crank-Nicolson")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("CN_comparacio.png", dpi=300)
plt.show()


# Per acabar, fem un resúm final

print("\n====================================")
print("RESUM DELS RESULTATS")
print("====================================")

for resultat in resultats:

    print(
        f'dt = {resultat["dt"]:.6f} | '
        f'Error maxim = {np.max(resultat["error"]):.8e}'
    )

print("\nSimulacions finalitzades correctament.")

