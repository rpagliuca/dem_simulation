#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time

# Simulation parameters taken from BIRWISCH et al., 2009

# Material
N = 100 # Number of grains
#RADIUS = 1.5E-4 # Radius of each grain --- (m)
RADIUS = 1.E-2 # Radius of each grain --- (m)
MASS= 7.63E-8 # Mass of each grain --- (kg)
MU = 0.5 # Drag coefficient --- (dimensionless)
MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
E_TILDE = 1.0E7 # Young's modulus / (1 - nu^2) --- (Pa)
GAMMA_R = 1.0E6 # Gamma / R --- (Pa . s/m)
KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
MU_W = 0.15 # --- (dimensionless)

# Dye geometry
H = 15.E-2 # dye height (in m)
L = 2.E-2 # dye length (in m)

# Physical and math constants
G = 9.81 # gravity
PI = np.pi

# Misc parameters
T0 = 0. # Initial time --- (s)
TF = 1. # End time --- (s)
STEPS = 200 # Number of steps --- (integer)
DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
DIMENSIONS = 2. # Number of degrees of freedoms (x,y => 2; x,y,z =>3)

# Derivative constants
GRAVITAL_FORCE = MASS*G
EFFECTIVE_RADIUS = (RADIUS*RADIUS)/(RADIUS+RADIUS)

# Constants to be used with the data matrices
X = 0
Y = 1
VX = 2
VY = 3
FX = 4
FY = 5

# Initializing matrices
EMPTY_MATRIX = np.zeros((N, DIMENSIONS*3.)) # Each matrix has 2 dimensions x 3 types of values (X, Y, VX, VY, FX, FY)
current_matrix = np.copy(EMPTY_MATRIX)
last_matrix = np.copy(EMPTY_MATRIX)
current_matrix[:,X:Y+1] = np.random.random_sample((N, DIMENSIONS))

def add_gravity_force(current_matrix):
    current_matrix[:,5] =- GRAVITAL_FORCE
    return current_matrix

def calculate_contact_forces(other_particle, current_particle):
    # Distance between particles
    contact_forces = np.zeros(2)
    distance = ((current_particle[X]-other_particle[X])**2. + (current_particle[Y]-other_particle[Y])**2.)**0.5
    # Fix division by zero
    if distance == 0:
        distance = 1e-20
    radial_unitary_vector = (other_particle[X:Y+1] - current_particle[X:Y+1])/distance
    # There is contact forces only if particles are closer than the sum of their radii
    if distance <= 2.*RADIUS:

        # Repulsion force
        contact_forces += 2./3. * E_TILDE * EFFECTIVE_RADIUS**0.5 * distance**(3./2.) * radial_unitary_vector

        # Cohesion force is null for the single sphere case, because w is 0
        # contact_forces += - (4 * PI * w * E_TILDE)**0.5 * EFFECTIVE_RADIUS**(3./4.) * distance**(3./4.) * radial_unitary_vector

        # Damping force
        contact_forces += - GAMMA_R * RADIUS * (EFFECTIVE_RADIUS * distance)**0.5 * np.dot( (other_particle[VX:VY+1] - current_particle[VX:VY+1]), radial_unitary_vector ) * radial_unitary_vector

    return contact_forces

def region_of_contact_forces(current_matrix, current_particle_index, max_neighbour_index):
    # Particles i through max_neighbour_index are the ones which may interact (based solely on x distance)
    region_of_interest = current_matrix[current_particle_index+1:max_neighbour_index+1, :]
    if region_of_interest.size > 0:
        contact_forces = np.apply_along_axis(calculate_contact_forces, Y, region_of_interest, current_matrix[current_particle_index,:])
        region_of_interest[:,FX:FY+1] += contact_forces 
        current_matrix[current_particle_index, FX:FY+1] += -contact_forces.sum(axis=0)
    return current_matrix

def add_contact_forces(current_matrix):
    # Now we iterate over every particle, only accounting other particles which x_i - x_j <= RADIUS
    # Last particle shouldn't interact with any other. It already interacted with the previous ones.
    for i in range (0, N-1):
        # np.argmax returns the minimum index wich satisfies some arbitrary condition
        # TODO: improve the following line, because argmax does not stop at first ocurrence
        max_neighbour_offset = np.argmax(current_matrix[i+1:, X] > current_matrix[i, X] + RADIUS)
        max_neighbour_index = i + max_neighbour_offset
        if max_neighbour_offset == 0:
            max_neighbour_index = current_matrix[:, X].size-1
        current_matrix = region_of_contact_forces(current_matrix, i, max_neighbour_index)
    return current_matrix

def add_stoke_air_drag_force(current_matrix):
    return current_matrix
    return forces_matrix - 6 * PI * MU_A * RADIUS * velocities_matrix

def calculate_velocities(current_matrix, last_matrix):
    current_matrix[:,VX:VY+1] = last_matrix[:,VX:VY+1] + (current_matrix[:,FX:FY+1] + last_matrix[:,FX:FY+1])/(2*MASS) * DT
    return current_matrix

def calculate_positions(current_matrix, last_matrix):
    current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + last_matrix[:,FX:FY+1]/(2*MASS) * DT**2
    # Boundary conditions
    current_matrix[:,X] = np.clip(current_matrix[:,X], 0, L)
    current_matrix[:,Y] = np.clip(current_matrix[:,Y], 0, H)
    return current_matrix

# Initializing matplotlib scatter plot data
scatterPoints = None
plt.ion()
plt.show()
plt.axis([0, L, 0, H])

# Solving steps
for step in range(1, STEPS+1):

    # First, we sort by x position to optimize contact forces
    current_matrix = current_matrix[current_matrix[:,X].argsort()]

    # Store current matrix as last matrix
    last_matrix = np.copy(current_matrix)

    # Reset matrix of forces
    current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

    # Apply gravity
    current_matrix = add_gravity_force(current_matrix)
    current_matrix = add_stoke_air_drag_force(current_matrix)
    current_matrix = add_contact_forces(current_matrix)

    # Calculate velocities from force
    current_matrix = calculate_velocities(current_matrix, last_matrix)

    # Calculate position from force and velocity
    current_matrix = calculate_positions(current_matrix, last_matrix)

    # Draw animation frame
    if scatterPoints:
        scatterPoints.remove()
    scatterPoints = plt.scatter(current_matrix[:, X], current_matrix[:,Y], s=100)
    plt.draw()
    plt.pause(1.0E-10)

print 'Finished.'
