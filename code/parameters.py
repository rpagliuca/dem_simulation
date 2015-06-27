#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit

# Reproduce random results for debugging
np.random.seed(1)

# Flags
realtimePlot = True
stepPlotFlag = False
lastPlotFlag = False

# Physical and math constants
G = 9.81 # gravity
PI = np.pi
INFINITY = 1.0E50

# Shoe dimensions
SH = 5.E-2 # shoe height (in m)
# Multiplicator for random generation
SH_MULTIPLICATOR = 2.
#L = 4.E-2 # dye length (in m)
SL = SH # show lengths (in m)

# Dye dimensions
DH = 2.E-2
DL = DH

# Particle size
#RADIUS = 5.E-4 # Radius of each grain --- (m)
RADIUS = 1.5E-3 # Radius of each grain --- (m)
scatterPlotPointSize = 1.0E8 * RADIUS**2

# Number of particles needed to represent the dye (they overlap a little bit ~ 1.6 instead of 2)
NUMBER_PARTICLES_BOTTOM_WALL = np.ceil(SL/(RADIUS*1.6))
NUMBER_PARTICLES_SIDE_WALL = np.ceil(SH/(RADIUS*1.6))
NUMBER_PARTICLES_DYE_BOTTOM_WALL = np.ceil(DL/(RADIUS*1.6))
NUMBER_PARTICLES_DYE_SIDE_WALL = np.ceil(DH/(RADIUS*1.6))

# Simulation parameters taken from BIRWISCH et al., 2009
# Material
N = 100
N += int(2*NUMBER_PARTICLES_BOTTOM_WALL + 2*NUMBER_PARTICLES_SIDE_WALL + NUMBER_PARTICLES_DYE_BOTTOM_WALL + 2*NUMBER_PARTICLES_DYE_SIDE_WALL) # Number of grains
MASS= 7.63E-8 # Mass of each grain --- (kg)
MU = 0.5 # Drag coefficient --- (dimensionless)
MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
MU_W = 0.15 # --- (dimensionless)

# Repulsion parameters
#E_TILDE = 1.0E-30 # Young's modulus / (1 - nu^2) --- (Pa)
E_TILDE = 1.0E32 * 1.0E-30 # Young's modulus / (1 - nu^2) --- (Pa)
#E_TILDE = 0.

# Viscous parameter
GAMMA_R = 1.0E-4 * 1.0E6 # Gamma / R --- (Pa . s/m)

# Friction parameter
GBPM_GAMMA = 1.E-6
#GBPM_GAMMA = 0

# Material of wall
WALL_MASS = INFINITY

# Misc parameters
T0 = 0. # Initial time --- (s)
#TF = 3. # End time --- (s)
STEPS = 2000 # Number of steps --- (integer)
#STEPS = 5 # Number of steps --- (integer)
#DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
DT = 5.E-4
TF = T0 + STEPS*DT
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
WT = 8 # Detailed type for walls (0 -> movable, 1 -> fixed)
