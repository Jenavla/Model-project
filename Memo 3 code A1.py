import numpy as np
import matplotlib.pyplot as plt

# Parameters
m = 4.1e-6  # massa
gamma = 7.3603e-7  # dempingcoëfficiënt
k = 0.1619  # veerconstante
Fmax = 60e-9  # amplitude van de oscillerende kracht op resonantie
dt = 0.0001  # tijdstapgrootte

# Bereken resonantiefrequentie
resonantiefrequentie = np.sqrt(k / m) / (2 * np.pi)
print('Resonantiefrequentie:', resonantiefrequentie)

# Tijdstappen
tmax = 100 * (1 / resonantiefrequentie)
num_steps = int(tmax / dt)
t = np.linspace(0, tmax, num_steps)

# Krachten
F = Fmax * np.sin(2 * np.pi * resonantiefrequentie * t) 

# Arrays voor positie en snelheid
x = np.zeros(num_steps)
v = np.zeros(num_steps)

# startposities
x[0] = 0.0
v[0] = 0.0

# Numerieke oplossing van de differentiaalvergelijking
for i in range(1, num_steps):
    x[i] = x[i - 1] + dt * v[i - 1]
    v[i] = v[i - 1] + dt * (F[i - 1] - gamma * v[i - 1] - k * x[i - 1] - 0.1 * v[i - 1]**3) / m 

# Berekenen amplitude bij statische kracht
x_Fmax = Fmax / k
print('Amplitude bij statische kracht:', x_Fmax)

# Bepalen van de maximale amplitude na stabilisatie
stabilization_steps = int(0.8 * num_steps)  # Aantal stappen na stabilisatie
max_amplitude = np.max(np.abs(x[stabilization_steps:]))
print('Amplitude na stabilisatie:', max_amplitude)

# Parameters voor tuning-kromme
num_frequencies = 100
frequency_range = 0.1 * resonantiefrequentie  # Range rond resonantiefrequentie om te onderzoeken

# Array van frequenties voor tuning-kromme
frequencies = np.linspace(resonantiefrequentie - frequency_range, resonantiefrequentie + frequency_range, num_frequencies)

# Array voor amplitudes voor tuning-kromme
amplitudes = np.zeros(num_frequencies)

# Berekenen van amplitudes bij verschillende frequenties voor tuning-kromme
for i in range(num_frequencies):
    F = Fmax * np.sin(2 * np.pi * frequencies[i] * t)

    z = np.zeros(num_steps)
    v = np.zeros(num_steps)
    z[0] = 0.0
    v[0] = 0.0

    for j in range(1, num_steps):
        z[j] = z[j-1] + dt * v[j-1]
        v[j] = v[j-1] + dt * (F[j-1] - gamma * v[j-1] - k * z[j-1]) / m

    stabilization_steps = int(0.8 * num_steps)
    max_amplitude = np.max(np.abs(z[stabilization_steps:]))
    amplitudes[i] = max_amplitude

# Bepaal FWHM
half_max_amplitude = max_amplitude / 2
half_max_indices = np.where(amplitudes >= half_max_amplitude)[0]
fwhm = frequencies[half_max_indices[-1]] - frequencies[half_max_indices[0]]
print('FWHM:', fwhm)

# Bepalen van de Q-factor
versterkingsfactor = resonantiefrequentie / fwhm
print('Versterkingsfactor (Q-factor):', versterkingsfactor)

# factor V0/FWHM
factor_fwhm = resonantiefrequentie / fwhm
print('Factor 𝜈₀/FWHM:', factor_fwhm)

# Plot tuning curve
fig1, ax1 = plt.subplots()
ax1.plot(frequencies, amplitudes)
ax1.set_xlabel('Frequentie(Hz)')
ax1.set_ylabel('Amplitude(m)')
ax1.set_title('Tuning-kromme van de drive mode')
ax1.grid(True)

# Bereken FWHM
fwhm = frequencies[half_max_indices[-1]] - frequencies[half_max_indices[0]]

# Voeg tekstlabel toe voor FWHM
ax1.annotate('FWHM: {:.2f}'.format(fwhm), xy=(frequencies[half_max_indices[0]], half_max_amplitude),
             xytext=(frequencies[half_max_indices[0]] + 0.2, half_max_amplitude*1.2))

plt.show()

# Plot uitwijking
fig2, ax2 = plt.subplots()
ax2.plot(t, x)
ax2.set_xlabel('t(s)')
ax2.set_ylabel('x(t)')
ax2.set_title('Uitwijking van de massa')
ax2.grid(True)

plt.show()