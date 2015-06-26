import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit
from parameters import *

# Initializing matrices
EMPTY_MATRIX = np.zeros((N, DIMENSIONS*3. + 3.)) # Each matrix has 2 dimensions x 3 types of values (X, Y, VX, VY, FX, FY) plus additional columns: Mass (M), Type (T) and Wall type (WT)
current_matrix = np.copy(EMPTY_MATRIX)
last_matrix = np.copy(EMPTY_MATRIX)
current_matrix[:,M] = MASS
current_matrix[:,T] = 1
current_matrix[:,WT] = 0 # not a wall

# Particles positions
current_matrix[:,X] = 2.0*RADIUS + np.random.random_sample(N)*(SL - 4.0*RADIUS)
current_matrix[:,Y] = 2.0*RADIUS + np.random.random_sample(N)*(SH - 4.0*RADIUS) * SH_MULTIPLICATOR

# Now we'll turn some of the particles into wall particles
# Shoe Bottom wall
offset = 0
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, Y] = 0
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, X] = np.linspace(0, SL, num=NUMBER_PARTICLES_BOTTOM_WALL)
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, M] = WALL_MASS
current_matrix[offset:offset+NUMBER_PARTICLES_BOTTOM_WALL, WT] = 1 # static wall

# Shoe Left wall
offset = offset+NUMBER_PARTICLES_BOTTOM_WALL
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, X] = 0
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, Y] = np.linspace(0, SH, num=NUMBER_PARTICLES_SIDE_WALLS)
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, M] = WALL_MASS
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, WT] = 2 # movable wall

# Shoe Right wall
offset = offset+NUMBER_PARTICLES_SIDE_WALLS
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, X] = SL
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, Y] = np.linspace(0, SH, num=NUMBER_PARTICLES_SIDE_WALLS)
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, M] = WALL_MASS
current_matrix[offset:offset+NUMBER_PARTICLES_SIDE_WALLS, WT] = 2 # movable wall

# Shoe Bottom wall
offset = offset+NUMBER_PARTICLES_SIDE_WALLS
current_matrix[offset:offset+NUMBER_PARTICLES_DYE_BOTTOM_WALL, Y] = -DH
current_matrix[offset:offset+NUMBER_PARTICLES_DYE_BOTTOM_WALL, X] = SL + np.linspace(0, DL, num=NUMBER_PARTICLES_DYE_BOTTOM_WALL)
current_matrix[offset:offset+NUMBER_PARTICLES_DYE_BOTTOM_WALL, M] = WALL_MASS
current_matrix[offset:offset+NUMBER_PARTICLES_DYE_BOTTOM_WALL, WT] = 1 # movable wall

# Flag particles as type 0 (wall particle)
offset = offset+NUMBER_PARTICLES_DYE_BOTTOM_WALL
current_matrix[0:offset, T] = 0
