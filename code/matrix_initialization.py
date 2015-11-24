# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
import init_overlap_fix
import dye_profiles.stepped as dye_profile
import parameters as p # Load all global variables from parameters

def matrix_initialization():

    print 'Beggining matrix initialization...'

    dye_profile.init()
    p.N_WALL = dye_profile.calculateN() # Number of particles needed for the walls
    print p.N_WALL
    p.N = p.DESIRED_N_PARTICLES + p.N_WALL # Number of grains for initial matrix

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

    dye_profile.draw(current_matrix)

    print 'Finished matrix initialization...'

    current_matrix = init_overlap_fix.init_overlap_fix(current_matrix)

    p.current_matrix = current_matrix
