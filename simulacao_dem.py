#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit

# Reproduce random results for debugging
np.random.seed(6)

# Flags
realtimePlot = False
log = False
stepPlotFlag = False

if log:
    print 'Starting...'

# Physical and math constants
G = 9.81 # gravity
PI = np.pi
INFINITY = 1.0E50

# Dye dimensions
H = 10.E-2 # dye height (in m)
# Multiplicator for random generation
H_MULTIPLICATOR = 8
#L = 4.E-2 # dye length (in m)
L = H # dye length (in m)

# Particle size
RADIUS = 5.E-3 # Radius of each grain --- (m)
#RADIUS = 1.5E-4 # Radius of each grain --- (m)

# Number of particles needed to represent the dye (they overlap a little bit ~ 1.6 instead of 2)
NUMBER_PARTICLES_BOTTOM_WALL = np.ceil(L/(RADIUS*1.6))
NUMBER_PARTICLES_SIDE_WALLS = np.ceil(H/(RADIUS*1.6))

# Simulation parameters taken from BIRWISCH et al., 2009
# Material
N = 200 
N += int(NUMBER_PARTICLES_BOTTOM_WALL + NUMBER_PARTICLES_SIDE_WALLS*2) # Number of grains
MASS= 7.63E-8 # Mass of each grain --- (kg)
MU = 0.5 # Drag coefficient --- (dimensionless)
MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
E_TILDE = 1.0E-30 # Young's modulus / (1 - nu^2) --- (Pa)
GAMMA_R = 1.0E6 # Gamma / R --- (Pa . s/m)
KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
MU_W = 0.15 # --- (dimensionless)

# Material of wall
WALL_MASS = INFINITY

# Misc parameters
T0 = 0. # Initial time --- (s)
TF = 1. # End time --- (s)
STEPS = 10000 # Number of steps --- (integer)
#STEPS = 5 # Number of steps --- (integer)
DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
if log:
    print "DT: ", DT
DIMENSIONS = 2. # Number of degrees of freedoms (x,y => 2; x,y,z =>3)

# Derivative constants
EFFECTIVE_RADIUS = (RADIUS*RADIUS)/(RADIUS+RADIUS)

# Constants to be used as columns on the data matrices
X = 0 # Positon on x
Y = 1 # Position on y
VX = 2 # Velocity x-component
VY = 3 # Velocity y-component
FX = 4 # Force x-component
FY = 5 # Force y-component
M = 6 # Mass
T = 7 # Type (particle -> 1, wall -> 0)

if log:
    print 'Initializing matrices...'

# Initializing matrices
EMPTY_MATRIX = np.zeros((N, DIMENSIONS*3. + 2.)) # Each matrix has 2 dimensions x 3 types of values (X, Y, VX, VY, FX, FY) plus additional columns: Mass
current_matrix = np.copy(EMPTY_MATRIX)
last_matrix = np.copy(EMPTY_MATRIX)
current_matrix[:,M] = MASS
current_matrix[:,T] = 1

# Particles positions
current_matrix[:,X] = np.random.random_sample(N)*L
current_matrix[:,Y] = np.random.random_sample(N)*H*H_MULTIPLICATOR

# Now we'll turn some of the particles into wall particles
# Bottom wall
offset = 0
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, Y] = 0
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, X] = np.linspace(0, L, num=NUMBER_PARTICLES_BOTTOM_WALL)
if log:
    print NUMBER_PARTICLES_BOTTOM_WALL
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, M] = WALL_MASS

# Left wall
offset = offset+NUMBER_PARTICLES_BOTTOM_WALL
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, X] = 0
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, Y] = np.linspace(0, H, num=NUMBER_PARTICLES_SIDE_WALLS)
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, M] = WALL_MASS

# Right wall
offset = offset+NUMBER_PARTICLES_SIDE_WALLS
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, X] = L
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, Y] = np.linspace(0, H, num=NUMBER_PARTICLES_SIDE_WALLS)
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, M] = WALL_MASS

# Flag particles as type 0 (wall particle)
offset = offset+NUMBER_PARTICLES_SIDE_WALLS
current_matrix[0:offset, T] = 0

if log:
    print 'Done!'

# Dye wall

# Source: http://stackoverflow.com/a/29799815/1501575
# Pre-compiled function to find first element of array greater than
@jit(nopython=True)
def find_first_item_greater_than(vec, value):
    for i in xrange(len(vec)):
        if vec[i] > value:
            return i
    return -1

def apply_forces(current_matrix):
    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= 2*RADIUS
    # Last particle shouldn't interact with any other. It has already interacted with the previous ones.
    for i in range (0, N-1):

        #print '---'
        #print 'i: ', i
        # np.argmax returns the minimum index wich satisfies some arbitrary condition, but it is way slower than using numba pre-compiled functions
        #max_neighbour_offset = np.argmax(current_matrix[i+1:, Y] > current_matrix[i, Y] + 2*RADIUS)
        max_neighbour_offset = find_first_item_greater_than(current_matrix[i+1:, Y], current_matrix[i, Y] + 2*RADIUS)
        max_neighbour_index = i+1 + max_neighbour_offset

        if max_neighbour_offset == -1:
            # Skip loop,there is no possible interaction for tihs particle
            continue

        # Particles i through max_neighbour_index are the ones which may interact (based solely on Y distance)
        row_of_interest = current_matrix[i+1:max_neighbour_index+1, :]

        # Now we remove those particles which Y position is greater than radius
        possible_interactions = np.less_equal(row_of_interest[:,X], current_matrix[i, X] + 2*RADIUS)
        cell_of_interest = row_of_interest[possible_interactions, :]

        #contact_forces = 0

        # If there is any possible neighbour
        if cell_of_interest.size > 0:

            # Get distance of possible neighbours
            distances = np.sqrt(np.square(cell_of_interest[:,X] - current_matrix[i,X]) + np.square(cell_of_interest[:,Y] - current_matrix[i, Y]))

            # Generate unitary vector
            radial_unitary_vector = ((cell_of_interest[:, X:Y+1] - current_matrix[i, X:Y+1]).transpose() / (distances + 1.E-20)).transpose()

            # Discard distances greater than 2*RADIUS
            deformations = (2*RADIUS - distances).clip(min=0)

            # Add contact forces

            # Force 1 => Repulsion force
            #contact_forces = (2./3. * E_TILDE * EFFECTIVE_RADIUS**0.5 * deformations**(3./2.) * radial_unitary_vector.transpose()).transpose()
            contact_forces = (2*7.48E-7 * 5.0E3 * deformations * radial_unitary_vector.transpose()).transpose()
            #contact_forces = (2*7.48E-7 * 5.0E0 / (distances + 1.E-20)**2 * deformations * radial_unitary_vector.transpose()).transpose()

            # Gravity applies force of approx. 9.81*7.63E-8 => 7.48E-7 NEWTON
            # Repulsion should apply more or less the same ammount (maybe double)

            # Force 2 => Attraction (cohesion)
            # Non-existent for single spheres

            # Force 3 => Viscous dissipation
            #contact_forces += - (GAMMA_R * RADIUS * np.sqrt(EFFECTIVE_RADIUS * deformations).transpose() * np.einsum( 'ij, ij->i', (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()
            contact_forces += - 1E-5 * (GAMMA_R * RADIUS * np.sqrt(EFFECTIVE_RADIUS * deformations).transpose() * np.einsum( 'ij, ij->i', (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()

            # Force 4 => Friction force
            # To be modelled
            # I'm not modelling as Cundall, but as Haff and Werner
            relative_velocities = (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1])
            tangent_relative_velocities = relative_velocities - (relative_velocities.transpose() - np.einsum('ij, ij->i', relative_velocities, radial_unitary_vector)).transpose() * radial_unitary_vector
            tangent_relative_velocities_module = np.sqrt(tangent_relative_velocities[:, 0]**2 + tangent_relative_velocities[:, 1]**2) + 1.E-20
            tangent_unitary_vector = (tangent_relative_velocities.transpose() / tangent_relative_velocities_module).transpose()
            GBPM_GAMMA = 1.E-6
            contact_forces += - GBPM_GAMMA * (tangent_relative_velocities.transpose() * deformations).transpose()

            # Write contact_forces to row_of_interest view based on possible_interactions items (indeces)
            row_of_interest[possible_interactions,FX:FY+1] += contact_forces 

            # Add forces and apply its negative sum to current particle (Newton's Second Law of motion)
            current_matrix[i, FX:FY+1] += -np.einsum('ij->j', contact_forces) 
            
    # Apply gravity and Stoke's air drag (except for wall particles by multiplying by column T)
    # Force 5 => Gravity
    # Force 6 => Stoke's air drag
    current_matrix[:, FX:FY+1] += ((-6 * PI * MU_A * RADIUS * current_matrix[:, VX:VY+1] + (np.outer(current_matrix[:, M]*G, np.array((0, -1))))).transpose() * current_matrix[:, T]).transpose()

    return current_matrix

# Initializing matplotlib scatter plot data
stepPlot = False
lastPlot = False
if realtimePlot:
    scatterPoints = None
    plt.ion()
    plt.show()
    plt.axis([0 - L*0.05, L + L*0.05, 0 - 2*H*0.05, 2*H + 2*H*0.05])
    plt.axes().set_aspect('equal')

if log:
    print 'Initializing steps...'

# Author: dr jimbob
# Source: http://stackoverflow.com/a/5009578/1501575
def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = pylab.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

# Solving steps
for step in range(1, STEPS+1):
    if log:
        print 'Beggining step ', str(step), '...'

    # First, we sort by y position to optimize contact forces
    if log:
        print 'Sorting matrix...'

    current_matrix = current_matrix[current_matrix[:,Y].argsort()]

    if log:
        print 'Done sorting matrix.'

    # Store current matrix as last matrix
    last_matrix = np.copy(current_matrix)

    # Reset matrix of forces
    current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

    # Apply forces
    if log:
        print 'Applying forces...'

    current_matrix = apply_forces(current_matrix)

    if log:
        print 'Done applying forces.'


    # Calculate velocities from force
    current_matrix[:,VX:VY+1] = last_matrix[:,VX:VY+1] + ((current_matrix[:,FX:FY+1] + last_matrix[:,FX:FY+1]).transpose()*current_matrix[:,T]/(2*current_matrix[:,M])).transpose() * DT

    #print current_matrix[:,VX:T+1]

    # Calculate position from force and velocity
    #current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + (last_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2
    current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + current_matrix[:,VX:VY+1]*DT + (current_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2

    # Draw animation frame
    if step == STEPS:
        lastPlot = True
    if stepPlotFlag and (step == STEPS or step == np.round(STEPS/2) or step == np.round(STEPS/4)):
        stepPlot = True
    if not realtimePlot and (stepPlot or lastPlot):
        scatterPoints = None
        plt.ion()
        plt.show()
        plt.axis([0 - L*0.05, L + L*0.05, 0 - 2*H*0.05, 2*H + 2*H*0.05])
        plt.axes().set_aspect('equal')
        realtimePlot = True

    if realtimePlot:
        if scatterPoints:
            scatterPoints.remove()
        scatterPoints = plt.scatter(current_matrix[:, X], current_matrix[:,Y], s=900, facecolors='none')
        if stepPlot or lastPlot:
            plt.ioff()
            plt.show()
        else:
            plt.draw()
            plt.pause(1.0E-10)

    if stepPlot:
        realtimePlot = False
        stepPlot = False

    print 'Done step ', str(step), '.' 

if log:
    print current_matrix
    print 'Finished.'
