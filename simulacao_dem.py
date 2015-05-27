#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time

# Simulation parameters taken from BIRWISCH et al., 2009

# Material
N = 20. # Number of grains
RADIUS = 1.5E-4 # Radius of each grain --- (m)
MASS= 7.63E-8 # Mass of each grain --- (kg)
MU = 0.5 # Drag coefficient --- (dimensionless)
MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
E_TILDE = 1.0E7 # Young's modulus / (1 - nu^2) --- (Pa)
GAMMA_R = 1.0E6 # Gamma / R --- (Pa . s/m)
KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
MU_W = 0.15 # --- (dimensionless)

# Dye geometry
H = 1. # dye height
L = 1. # dye length

# Physical and math constants
G = 9.81 # gravity
PI = np.pi

# Misc parameters
T0 = 0. # Initial time --- (s)
TF = 1 # End time --- (s)
STEPS = 20 # Number of steps --- (integer)
DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
DIMENSIONS = 2. # Number of degrees of freedoms (x,y => 2; x,y,z =>3)

# Derivative constants
GRAVITAL_FORCE = MASS*G
EFFECTIVE_RADIUS = (RADIUS*RADIUS)/(RADIUS+RADIUS)

# Initializing matrices
EMPTY_MATRIX = np.zeros((N, DIMENSIONS))
forces_matrix = np.copy(EMPTY_MATRIX)
last_forces_matrix = np.copy(EMPTY_MATRIX)
last_positions_matrix = np.copy(EMPTY_MATRIX)
last_velocities_matrix = np.copy(EMPTY_MATRIX)
velocities_matrix = np.copy(EMPTY_MATRIX)
positions_matrix = np.random.random_sample((N, DIMENSIONS))
print positions_matrix

def add_gravity_force(forces_matrix):
    forces_matrix[:,1] -= GRAVITAL_FORCE
    return forces_matrix

def add_repulsion_force(forces_matrix, positions_matrix):
    # First, we sort by x position so we can easily ignore lots of interactions
    return forces_matrix

def add_stoke_air_drag_force(forces_matrix, velocities_matrix):
    return forces_matrix - 6 * PI * MU_A * RADIUS * velocities_matrix

def calculate_velocities(forces_matrix, last_forces_matrix, last_velocities_matrix):
    return last_velocities_matrix + (forces_matrix + last_forces_matrix)/(2*MASS) * DT

def calculate_positions(last_forces_matrix, last_positions_matrix, last_velocities_matrix):
    positions_matrix = last_positions_matrix + last_velocities_matrix*DT + last_forces_matrix/(2*MASS) * DT**2
    # Boundary conditions
    positions_matrix[:,0] = np.clip(positions_matrix[:,0], 0, L)
    positions_matrix[:,1] = np.clip(positions_matrix[:,1], 0, H)
    return positions_matrix

plt.ion()
plt.show()
plt.axis([0, L, 0, H])

# Initializing matplotlib scatter plot data
scatterPoints = None

# Solving steps
for step in range(1, STEPS+1):

    # Reset matrix of forces
    last_forces_matrix = np.copy(forces_matrix)
    forces_matrix = np.copy(EMPTY_MATRIX)

    # Apply gravity
    forces_matrix = add_gravity_force(forces_matrix)
    forces_matrix = add_stoke_air_drag_force(forces_matrix, velocities_matrix)

    # Calculate velocities from force
    last_velocities_matrix = np.copy(velocities_matrix)
    velocities_matrix = calculate_velocities(forces_matrix, last_forces_matrix, last_velocities_matrix)

    # Calculate position from force and velocity
    last_positions_matrix = np.copy(positions_matrix)
    positions_matrix = calculate_positions(last_forces_matrix, last_positions_matrix, last_velocities_matrix)

    # Draw animation frame
    if scatterPoints:
        scatterPoints.remove()
    scatterPoints = plt.scatter(positions_matrix[:,0], positions_matrix[:,1], s=500)
    plt.draw()
    plt.pause(1.0E-10)
