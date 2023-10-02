import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Week 1 code --------------------------------------------------------------------------------------------------------------------

# lees de data from file
with open('posities_1_Team_A1.txt') as f:
    lines = f.readlines()
with open('posities_2_Team_A1.txt') as g:
    lines_2 = g.readlines()


# Maak twee kolommen van de text file
data = []
for line in lines:
    values = line.strip().split()  # Split de data op basis van spaties en maakt van alles een float
    try:
        value_1 = float(values[0])  
        value_2 = float(values[1])  
        data.append([value_1, value_2])  
    except ValueError:
        print("Error: Could not convert string to float on dataset")
data_2 = []
for line_2 in lines_2:
    values = line_2.strip().split()  # Split de data op basis van spaties en maakt van alles een float
    try:
        value_3 = float(values[0])  
        value_4 = float(values[1])  
        data_2.append([value_3, value_4])  
    except ValueError:
        print("Error: Could not convert string to float on dataset")

# Maak de datalijst
df = pd.DataFrame(data, columns=['column 1', 'column 2'])
df_2 = pd.DataFrame(data_2, columns=['column 1', 'column 2'])

# Bereken de versnelling
tijd = df['column 1']
positie = df['column 2']
snelheid = np.gradient(positie, tijd)
versnelling = np.gradient(np.gradient(positie, tijd), tijd)

tijd_2 = df_2['column 1']
positie_2 = df_2['column 2']
snelheid_2 = np.gradient(positie_2, tijd_2)
versnelling_2 = np.gradient(np.gradient(positie_2, tijd_2), tijd_2)


# Week 2 code --------------------------------------------------------------------------------------------------------------

# Constanten
m = 1.3837e-6 # Massa
k = 30 # Veerconstante
kritische_demping = 4*m*k
print('kritische demping:', kritische_demping)
gamma = 0.032215 #dempingscoëfficient

# Stapgrootte
dt = tijd[1]
t = tijd
dT = tijd_2[1]
T = tijd_2

# Startpunten
x0 = 0   # Positie
v0 = 0  # Snelheid
x02 = 0   # Positie
v02 = 0  # Snelheid

# Lijsten voor de data
x = [x0]
v = [v0]
x_2 = [x02]
v_2 = [v02]

# Numeriek oplossen
for i in range(0, len(t)-1):
    x.append(x[-1] + v[-1]*dt)
    v.append(v[-1] + (versnelling[i]-k*x[-2]/m - gamma*v[-1]/m)*dt)
for i in range(0, len(T)-1):
    x_2.append(x_2[-1] + v_2[-1]*dT)
    v_2.append(v_2[-1] + (versnelling_2[i] -k*x_2[-2]/m - gamma*v_2[-1]/m)*dT)

# Bepaal de locatie van de pieken in de eerste dataset
peaks1 = np.where(np.diff(np.sign(np.diff(versnelling))) == -2)[0] + 1

# Bepaal de locatie van de pieken in de tweede dataset
peaks2 = np.where(np.diff(np.sign(np.diff(x))) == -2)[0] + 1

piekverschil = []

# Vergelijk de pieken tussen beide datasets
for peak1 in peaks1:
    closest_peak2 = np.argmin(np.abs(peaks2 - peak1))
    peak_difference = np.abs(tijd[peak1] - tijd_2[peaks2[closest_peak2]])
    piekverschil.append(peak_difference)

gem_verschil = sum(piekverschil)/len(piekverschil)

print("Het gemiddelde verschil tussen de pieken in dataset 1 en dataset 2 is:", gem_verschil)


# Plots
fig, ax0 = plt.subplots()
ax0.plot(t, x, label = 'Positie 1')
ax0.set_xlabel('Tijd(s)')
ax0.set_ylabel('Positie(m)')
ax1 = ax0.twinx()
ax1.plot(tijd, versnelling, label = 'Versnelling 1', color = 'orange')
ax1.set_ylabel('versnelling(m/s^2)')
fig.legend(labels=['Respons', 'Versnelling'])
ax0.set_title('Posities 1')

fig, ax2 = plt.subplots()
ax2.set_xlabel('Tijd(s)')
ax2.set_ylabel('Positie')
ax2.plot(T, x_2, label = 'Positie 2')
ax3 = ax2.twinx()
ax3.plot(tijd_2, versnelling_2, label = 'Versnelling 2', color = 'orange')
ax3.set_ylabel('Versnelling(m/s^2)')
ax2.set_title('Posities 2')
fig.legend(labels=['Respons', 'Versnelling'])

plt.show()