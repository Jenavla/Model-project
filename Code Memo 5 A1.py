import numpy as np
import matplotlib.pyplot as plt

# Parameters stap 1
d_drive = 2e-3  # Dikte drive
l_drive = 200e-6  # Lengte drive
w_drive = 3e-6  # Diepte drive
N_drive = 100  # Aantal windingen
V_drive = 1.5  # Spanning over de drive
V_max = 15  # Maximale spanning over de drive
A_drive = l_drive * w_drive  # Oppervlakte van een plaat
mu0 = 4 * np.pi * 1e-7  # Permeabiliteit van vacuüm
E0 = 8.854e-12  # Elektrische veldconstante
dt = 1e-7

# Parameters stap 2
m_drive = 4.3e-9 # massa drive
k_drive = 0.3820 # veerconstante drive
gamma_drive = 8.1053e-07 # Dempingscoëfficiënt drive

# Parameters stap 3
omega = 10  # Rotatiesnelheid

# Parameters stap 4
m_sense = 5.16e-9  # Massa sense
k_sense = 1.8334  # Veerconstante sense
gamma_sense = 1.945e-06  # Dempingscoëfficiënt sense

# Parameters stap 5
l_sense = 200e-9 # lengte sense
d_sense = 2e-9 # breedte sense
w_sense = 3e-9 # diepte sense
N_sense = 40 # aantal blauwe platen
Vdc = 15 # spanning door sense
A_sense = l_sense * w_sense # oppervlakte sense condensatoren

# --------------------

# Bereken resonantiefrequentie
resonantie_drive = np.sqrt(k_drive / m_drive) 
print('Resonantiefrequentie Drive mode:', resonantie_drive)

# Tijdstappen
tmax = 0.05
num_steps = int(tmax / dt)
t = np.linspace(0, tmax, num_steps)

# Fel berekenen
A = (N_drive + 2) * (E0 * w_drive / (2 * d_drive)) * (V_max**2 + 0.5 * V_drive**2) # constante kracht dus doet niet veel voor beweging
B = (N_drive + 2) * (E0 * w_drive / (d_drive)) * V_max * V_drive
C = (N_drive + 2) * (E0 * w_drive / (4 * d_drive)) * V_drive**2 # te klein

Fel = B * np.cos(resonantie_drive * t) 
print('gemiddelde elektrische kracht:', np.mean(Fel))

# Numerieke oplossing van de differentiaalvergelijking
x_drive = np.zeros(num_steps)

stap = t[1] - t[0]

a = (k_drive - 2 * m_drive / (stap ** 2)) / (m_drive / (stap ** 2) + gamma_drive / (2 * stap))
b = (m_drive / (stap ** 2) - gamma_drive / (2 * stap)) / (m_drive / (stap ** 2) + gamma_drive / (2 * stap))
c = Fel / (m_drive / (stap ** 2) + gamma_drive / (2 * stap))

x_drive[0] = 0.0
x_drive[1] = 0.0

for i in range(1, num_steps - 1):
    x_drive[i + 1] = c[i] - a * x_drive[i] - b * x_drive[i - 1]

# Berekenen corioliskracht
v_lijst = np.zeros(num_steps)
v_lijst[0] = 0.0

for i in range(0, len(x_drive)-1):
    v_lijst[i+1] = (x_drive[i+1]-x_drive[i])/(t[i+1]-t[i])
    
crossproduct = v_lijst*omega
F_coriolis = -2*m_drive*crossproduct
print('gemiddelde Coriolis kracht:', np.mean(F_coriolis))

# Numerieke oplossing van de differentiaalvergelijking x-sense
x_sense = np.zeros_like(t)

stap = t[1] - t[0]

d = (k_sense - 2 * m_sense / (stap ** 2)) / (m_sense / (stap ** 2) + gamma_sense / (2 * stap))
e = (m_sense / (stap ** 2) - gamma_sense / (2 * stap)) / (m_sense / (stap ** 2) + gamma_sense / (2 * stap))
f = F_coriolis / (m_sense / (stap ** 2) + gamma_sense / (2 * stap))

x_sense[0] = 0.0
x_sense[1] = 0.0

for i in range(1, len(t) - 1):
    x_sense[i + 1] = f[i] - d * x_sense[i] - e * x_sense[i - 1]

# Bereken ∆C
delta_C_lijst = []  # Lijst om de delta_C waarden op te slaan

for i in range(1, len(x_sense)):
    delta_d_sense = abs(x_sense[i] - x_sense[i - 1])

    delta_C = (E0 * (A_sense / d_sense)) * (delta_d_sense / d_sense)
    delta_C_lijst.append(delta_C)

print('Gemiddelde capaciteitsverandering:', sum(delta_C_lijst) / len(delta_C_lijst))

# Bereken de lading Q op de condensator op elk tijdstip
Q = (E0 * l_sense * w_sense * Vdc * N_sense) / (d_sense-x_sense)

# Bereken de stroom I door de afgeleide van Q te nemen
I_sense = np.zeros_like(Q)

for i in range(1, len(Q)-1):
    I_sense[i] = (Q[i+1]-Q[i])/num_steps

print("I_sense:", I_sense)

# Transmissiecoëfficiënten
N1 = np.max(x_drive) / np.max(Fel)
print('Transmissiecoëfficiënt 1:', N1, 'm/N')

N2 = np.max(F_coriolis) / np.max(x_drive)
print('Transmissiecoëfficiënt 2:', N2, 'N/m')

N3 = np.max(x_sense) / np.max(F_coriolis)
print('Transmissiecoëfficiënt 3:', N3, 'm/N')

N4 = np.max(delta_C_lijst) / np.max(x_sense)
print('Transmissiecoëfficiënt 4:', N4, 'F/m')

N5 = np.max(I_sense) / np.max(delta_C_lijst)
print('Transmissiecoëfficiënt 5:', N5, 'A/F')

# Plot x-drive
fig1, ax1 = plt.subplots()
ax1.plot(t, x_drive)
ax1.set_xlabel('t(s)')
ax1.set_ylabel('x(t)')
ax1.set_title('Verplaatsing van de x-drive')
ax1.grid(True)

# Plot Fel
fig2, ax2 = plt.subplots()
ax2.plot(t, Fel)
ax2.set_xlabel('t(s)')
ax2.set_ylabel('Fel(N)')
ax2.set_title('Elektrische kracht (Fel)')
ax2.grid(True)

# Plot F_coriolis
fig3, ax3 = plt.subplots()
ax3.plot(t, F_coriolis)
ax3.set_xlabel('t(s)')
ax3.set_ylabel('F_coriolis(N)')
ax3.set_title('Corioliskracht (F_coriolis)')
ax3.grid(True)

# Plot x-sense
fig4, ax4 = plt.subplots()
ax4.plot(t, x_sense)
ax4.set_xlabel('t(s)')
ax4.set_ylabel('x(t)')
ax4.set_title('Verplaatsing van de x-sense')
ax4.grid(True)

# Plot I_sense
fig5, ax5 = plt.subplots()
ax5.plot(t, I_sense)
ax5.set_xlabel('t(s)')
ax5.set_ylabel('I_sense(A)')
ax5.set_title('Stroom door de sense-condensator (I_sense)')
ax5.grid(True)

plt.show()