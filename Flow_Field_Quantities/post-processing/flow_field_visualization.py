#File:		flow_field_visualization.py
#Date:		19.07.2023
#Author:	Armand Sieber
#Tag:		Display the velocity and presssure fields


import numpy as np
import matplotlib.pyplot as plt

# ---- User-defined inputs ----
# -----------------------------

# Must be provided by user
file_name = 'free_surface_gamma08_0538.txt' # simulation output filename

# ---- Loads  solution vectors ----
# ---------------------------------
data = np.genfromtxt(file_name, delimiter=',', skip_header=1)

x_fluid = data[:,0] # x-coordinates
y_fluid = data[:,1] # y-coordinates
ux_fluid = data[:,3] # velocity x-direction --> not used here, but included for completeness
uy_fluid = data[:,4] # velocity y-direction --> not used here, but included for completeness
U_fluid = data[:,6] # velocity magnitude
P_fluid = data[:,7] # pressure

# Amount of grid point in x and y directions
Nx_domain = round(abs((np.max(x_fluid) - np.min(x_fluid))/(x_fluid[1] - x_fluid[0])) + 1)
Ny_domain = round(abs((np.max(y_fluid) - np.min(y_fluid))/(y_fluid[0]- y_fluid[Nx_domain])) + 1)

# Reshape solution vectors for compatibility with python "contourf" plots
x_domain = np.reshape(x_fluid, (Ny_domain, Nx_domain))
y_domain = np.reshape(y_fluid, (Ny_domain, Nx_domain))
velocity_field = np.reshape(U_fluid, (Ny_domain, Nx_domain))
pressure_field = np.reshape(P_fluid, (Ny_domain, Nx_domain))

# ---- Display flow fields ----
# -----------------------------

# Velocity contour plot
fig1 = plt.figure(1)
ax = fig1.add_subplot()
plt.axis('equal')
plt.title('Velocity field')
cont = ax.contourf(x_domain, y_domain, velocity_field, cmap ="turbo", levels = 100)
ax.contourf(-x_domain, y_domain, velocity_field, cmap ="turbo", levels = 100) # obtained by symmetry
cbar = plt.colorbar(cont)
cbar.set_label('U [-]')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim((-np.max(x_fluid), np.max(x_fluid)))

# Pressure contour plot
fig2 = plt.figure(2)
ax = fig2.add_subplot()
plt.axis('equal')
plt.title('Pressure field')
cont = ax.contourf(x_domain, y_domain, pressure_field, cmap ="turbo", levels = 100)
ax.contourf(-x_domain, y_domain, pressure_field, cmap ="turbo", levels = 100) # obtained by symmetry
cbar = plt.colorbar(cont)
cbar.set_label('P [-]')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim((-np.max(x_fluid), np.max(x_fluid)))

plt.show()
