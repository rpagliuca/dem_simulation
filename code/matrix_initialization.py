# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
import init_overlap_fix
import parameters as p # Load all global variables from parameters

def matrix_initialization():

    print 'Beggining matrix initialization...'

    # Initializing matrices
    EMPTY_MATRIX = np.zeros((p.N, p.DIMENSIONS*3. + 3.)) # Each matrix has 2 dimensions x 3 types of values (p.X, p.Y, Vp.X, VY, Fp.X, FY) plus additional columns: p.Mass (M), Type (T) and Wall type (p.WT)
    current_matrix = np.copy(EMPTY_MATRIX)
    last_matrix = np.copy(EMPTY_MATRIX)
    current_matrix[:,p.M] = p.MASS
    current_matrix[:,p.T] = 1
    current_matrix[:,p.WT] = 0 # not a wall

    # Particles positions
    current_matrix[:,p.X] = 2.0*p.RADIUS + np.random.random_sample(p.N)*(p.SL - 4.0*p.RADIUS)
    current_matrix[:,p.Y] = 2.0*p.RADIUS + np.random.random_sample(p.N)*(p.SH - 4.0*p.RADIUS) * p.SH_MULTIPLICATOR

    # Now we'll turn some of the particles into wall particles
    # Table Bottom wall 1
    offset = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.Y] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.X] = np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.WT] = 1 # static wall

    # Table Bottom Wall 2
    offset = offset+p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.Y] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.X] = p.SL + p.DL + np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.WT] = 1 # static wall

    # Shoe Left wall
    offset = offset+p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.X] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.WT] = 2 # movable wall

    # Shoe Right wall
    offset = offset+p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.X] = p.SL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.WT] = 2 # movable wall

    # Dye Bottom wall
    offset = offset+p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.Y] = -p.DH
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.X] = p.SL + np.linspace(0, p.DL, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.WT] = 1 # static

    # Dye Left Wall
    offset = offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.X] = p.SL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.WT] = 1 # static

    # Dye Right Wall
    offset = offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.X] = p.SL + p.DL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.WT] = 1 # static

    # Flag particles as type 0 (wall particle)
    offset = offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[0:offset, p.T] = 0

    print 'Finished matrix initialization...'

    current_matrix = init_overlap_fix.init_overlap_fix(current_matrix)

    p.current_matrix = current_matrix
