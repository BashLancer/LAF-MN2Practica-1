
import numpy as np
import matplotlib.pyplot as plt

# Paràmetres del problema
N = 101
L = 1.0
t_final = 0.025

dx = L / (N - 1)
x = np.linspace(0, L, N)

tol = 1e-12
max_iter = 100000

factors_dt = [1.0, 0.5]


# Solució analítica mitjançant la sèrie de Fourier
def solucio_analitica(x, t, termes=399):

    T = np.zeros_like(x, dtype=float)

    for n in range(1, termes + 1, 2):
        coef = 4 / (n * np.pi)

        T += (coef * (1 - np.exp(-n**2 * np.pi**2 * t))
              / (n**2 * np.pi**2) * np.sin(n * np.pi * x))

    return T


# Mètode de Gauss-Seidel per resoldre el sistema
def gauss_seidel(b, gamma, x0):

    sol = x0.copy()
    m = len(b)

    for k in range(max_iter):

        anterior = sol.copy()

        for j in range(m):

            if j == 0:
                esquerra = 0
            else:
                esquerra = sol[j - 1]

            if j == m - 1:
                dreta = 0
            else:
                dreta = sol[j + 1]

            sol[j] = (b[j] + gamma / 2 * (esquerra + dreta)) / (1 + gamma)

        error = np.max(np.abs(sol - anterior))

        if error < tol:
            return sol, k + 1

    raise RuntimeError("Gauss-Seidel no ha convergit")


# Esquema de Crank-Nicolson
def crank_nicolson(dt):

    gamma = dt / dx**2
    Nt = round(t_final / dt)

    assert np.isclose(Nt * dt, t_final)

    T = np.zeros(N)
    temps = np.zeros(Nt + 1)

    for n in range(Nt):

        # Terme independent del sistema
        b = (gamma / 2 * T[:-2]
             + (1 - gamma) * T[1:-1]
             + gamma / 2 * T[2:]
             + dt)

        # Fem servir el pas anterior com a valor inicial
        T0 = T[1:-1].copy()

        T_int, iteracions = gauss_seidel(b, gamma, T0)

        T_nova = np.zeros(N)
        T_nova[1:-1] = T_int

        T = T_nova
        temps[n + 1] = (n + 1) * dt

    return T, temps, gamma


# Càlcul per als dos passos temporals
resultats = []

T_exacta = solucio_analitica(x, t_final)

for factor in factors_dt:

    dt = factor * dx**2

    T_num, temps, gamma = crank_nicolson(dt)

    error = np.abs(T_num - T_exacta)

    error_max = np.max(error)
    error_mig = np.mean(error)

    print("\nCrank-Nicolson")
    print(f"Factor = {factor}")
    print(f"dx = {dx:.8f}")
    print(f"dt = {dt:.8f}")
    print(f"Gamma = {gamma:.4f}")
    print(f"Passos temporals = {len(temps) - 1}")
    print(f"Temperatura màxima = {np.max(T_num):.8f}")
    print(f"Error màxim = {error_max:.8e}")
    print(f"Error mitjà = {error_mig:.8e}")

    resultats.append({
        "factor": factor,
        "dt": dt,
        "T": T_num.copy(),
        "error": error.copy()
    })

    # Comparació amb la solució analítica
    plt.figure(figsize=(9, 5))

    plt.plot(x, T_num, label="Crank-Nicolson", linewidth=2)
    plt.plot(x, T_exacta, "--", label="Solució analítica", linewidth=2)

    plt.xlabel(r"$\hat{x}$")
    plt.ylabel(r"$\hat{T}$")
    plt.title(f"Crank-Nicolson: dt = {dt:.6f}")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(f"CN_temperatura_{factor}.png", dpi=300)
    plt.show()

    # Error absolut respecte a la solució analítica
    plt.figure(figsize=(9, 5))

    plt.plot(x, error, label="Error absolut", linewidth=2)

    plt.xlabel(r"$\hat{x}$")
    plt.ylabel("Error absolut")
    plt.title(f"Error de Crank-Nicolson: dt = {dt:.6f}")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(f"CN_error_{factor}.png", dpi=300)
    plt.show()


# Comparació dels dos passos temporals
plt.figure(figsize=(9, 5))

for resultat in resultats:

    plt.plot(x, resultat["T"], linewidth=2,
             label=f'dt = {resultat["dt"]:.6f}')

plt.xlabel(r"$\hat{x}$")
plt.ylabel(r"$\hat{T}$")
plt.title("Comparació dels dos casos de Crank-Nicolson")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("CN_comparacio.png", dpi=300)
plt.show()


# Resum dels errors
print("\nResum dels resultats")

for resultat in resultats:
    print(f'dt = {resultat["dt"]:.6f} | '
          f'Error màxim = {np.max(resultat["error"]):.8e}')

print("\nSimulacions finalitzades.")


