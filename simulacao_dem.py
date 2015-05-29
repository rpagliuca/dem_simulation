#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit

print 'Starting...'

# Reproduce random results for debugging
np.random.seed(1)

# Flags
plot = True

# Physical and math constants
G = 9.81 # gravity
PI = np.pi
INFINITY = np.inf

# Simulation parameters taken from BIRWISCH et al., 2009
# Material
N = 1000 # Number of grains
#RADIUS = 1.5E-4 # Radius of each grain --- (m)
RADIUS = 1.5E-4 # Radius of each grain --- (m)
MASS= 7.63E-8 # Mass of each grain --- (kg)
MU = 0.5 # Drag coefficient --- (dimensionless)
MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
E_TILDE = 10.0E7 # Young's modulus / (1 - nu^2) --- (Pa)
GAMMA_R = 1.0E6 # Gamma / R --- (Pa . s/m)
KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
MU_W = 0.15 # --- (dimensionless)

# Material of wall
WALL_MASS = INFINITY

# Dye geometry
H = 10.E-3 # dye height (in m)
L = 4.E-3 # dye length (in m)

# Misc parameters
T0 = 0. # Initial time --- (s)
TF = 0.5 # End time --- (s)
STEPS = 260000 # Number of steps --- (integer)
DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
print "DT: ", DT
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

print 'Initializing matrices...'

# Initializing matrices
EMPTY_MATRIX = np.zeros((N, DIMENSIONS*3.)) # Each matrix has 2 dimensions x 3 types of values (X, Y, VX, VY, FX, FY)
current_matrix = np.copy(EMPTY_MATRIX)
last_matrix = np.copy(EMPTY_MATRIX)
current_matrix[:,X:Y+1] = np.random.random_sample((N, DIMENSIONS))
current_matrix[:,X] = np.random.random_sample(N)*L
current_matrix[:,Y] = np.random.random_sample(N)*H

print 'Done!'

# Dye wall

# Source: http://stackoverflow.com/a/29799815/1501575
# Pre-compiled function to find first element of array greater than
@jit(nopython=True)
def find_first_item_greater_than(vec, item):
    for i in xrange(len(vec)):
        if item < vec[i]:
            return i
    return -1

def add_contact_forces(current_matrix):
    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= RADIUS
    # Last particle shouldn't interact with any other. It already interacted with the previous ones.
    for i in range (0, N-1):
        # np.argmax returns the minimum index wich satisfies some arbitrary condition, but it is way slower than using numba pre-compiled functions
        #max_neighbour_offset = np.argmax(current_matrix[i+1:, Y] > current_matrix[i, Y] + RADIUS)
        max_neighbour_offset = find_first_item_greater_than(current_matrix[i+1:, Y], current_matrix[i, Y] + RADIUS)
        max_neighbour_index = i + max_neighbour_offset
        if max_neighbour_offset == 0:
            max_neighbour_index = current_matrix[:, Y].size-1

        # Particles i through max_neighbour_index are the ones which may interact (based solely on X distance)
        region_of_interest = current_matrix[i+1:max_neighbour_index+1, :]
        # Now we remove those particles which Y position is greater than radius
        possible_interactions = np.less_equal(region_of_interest[:,X], current_matrix[i, X] + RADIUS)
        subregion_of_interest = region_of_interest[possible_interactions, :]

        # If there is any possible neighbour
        if subregion_of_interest.size > 0:

            # Get distance of possible neighbours
            distances = np.sqrt(np.square(subregion_of_interest[:,X] - current_matrix[i,X]) + np.square(subregion_of_interest[:,Y] - current_matrix[i, Y]))

            # Generate unitary vector
            radial_unitary_vector = ((subregion_of_interest[:, X:Y+1] - current_matrix[i, X:Y+1]).transpose() / (distances.transpose() + 1E-20)).transpose()

            # Discard distances greater than 2*RADIUS
            distances = distances * np.greater(2*RADIUS - distances, 0)

            # Add contact forces
            contact_forces = (2./3. * E_TILDE * EFFECTIVE_RADIUS**0.5 * distances**(3./2.) * radial_unitary_vector.transpose()).transpose()
            contact_forces += - (GAMMA_R * RADIUS * np.sqrt(EFFECTIVE_RADIUS * distances).transpose() * np.einsum( 'ij, ij->i', (subregion_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()
            region_of_interest[possible_interactions,FX:FY+1] += contact_forces 
            # Add forces and apply its negative value to current particle (Newton's Second Law of motion)
            current_matrix[i, FX:FY+1] += -np.einsum('ij->j', contact_forces)

    return current_matrix

# Initializing matplotlib scatter plot data
if plot:
    scatterPoints = None
    plt.ion()
    plt.show()
    plt.axis([0, L, 0, H])

print 'Initializing steps...'

# Solving steps
for step in range(1, STEPS+1):

    # First, we sort by x position to optimize contact forces
    current_matrix = current_matrix[current_matrix[:,X].argsort()]

    # Store current matrix as last matrix
    last_matrix = np.copy(current_matrix)

    # Reset matrix of forces
    current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

    # Apply gravity
    current_matrix[:,FY] =- GRAVITAL_FORCE
    # Apply stoke air drag 
    current_matrix[:, FX:FY+1] = current_matrix[:, FX:FY+1] - 6 * PI * MU_A * RADIUS * current_matrix[:, VX:VY+1]
    # Apply contact forces
    current_matrix = add_contact_forces(current_matrix)

    # Calculate velocities from force
    current_matrix[:,VX:VY+1] = last_matrix[:,VX:VY+1] + (current_matrix[:,FX:FY+1] + last_matrix[:,FX:FY+1])/(2*MASS) * DT

    # Calculate position from force and velocity
    current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + last_matrix[:,FX:FY+1]/(2*MASS) * DT**2
    # Boundary conditions
    current_matrix[:,X] = np.clip(current_matrix[:,X], 0, L)
    current_matrix[:,Y] = np.clip(current_matrix[:,Y], 0, H)

    # Draw animation frame
    if plot:
        if scatterPoints:
            scatterPoints.remove()
        scatterPoints = plt.scatter(current_matrix[:, X], current_matrix[:,Y], s=50)
        plt.draw()
        plt.pause(1.0E-10)
    
    print 'Step ', str(step), ' done.'

print current_matrix

print 'Finished.'
