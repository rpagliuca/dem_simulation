#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import time
import os

global np, realtimePlot, stepPlotFlag, lastPlotFlag, G, PI, INFINITY, load_simulation

# Reproduce random results for debugging
np.random.seed(1)

# Flags
load_simulation = False
realtimePlot = True
stepPlotFlag = False
lastPlotFlag = False

# Physical and math constants
G = 9.81 # gravity
PI = np.pi
INFINITY = 1.0E50

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
DIMENSIONS = 2. # Number of degrees of freedoms (x,y => 2; x,y,z =>3)

# If flag load_simulation is on
if load_simulation:

    import_path = "/home/rsantos/Desktop/simulacao_dem/output/simulation_RADIUS0.0008_DT0.0001_ETILDE800.0_GAMMAR200.0_GBPMGAMMA4e-05_N3208/step800"

    print "Loading state from:"
    print import_path
    print ""

    current_matrix = np.loadtxt(import_path + "/current_matrix.txt")
    SAVE_SESSION_STEP_INTERVAL, SAVE_SESSION_DIFFERENT_FILE_PER_STEP, SH, SL, SH_MULTIPLICATOR, DH, DL, N, RADIUS, scatterPlotPointSize, MASS, MU, MU_A, KAPPA_R, MU_W, shoe_velocity, step, T0, DT, STEPS, GBPM_GAMMA, GAMMA_R, E_TILDE = np.loadtxt(import_path + "/parameters.txt", unpack = True)

    print "N = " + str(N)

    SAVE_ENABLED = False

    # Convert some floats to int
    step = int(step)
    STEPS = int(STEPS)
    STEPS = 1000
    N = int(N)

    start_step = step + 1

else: # If flag load_simulation is off, then it is a new simulation and we should set the default parameters

    # Initial shoe velocity
    shoe_velocity = 0.0

    start_step = 1

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
    RADIUS = 2.E-4 # Radius of each grain --- (m)
    #RADIUS = 8.E-4 # Radius of each grain --- (m)
    scatterPlotPointSize = 1.0E8 * RADIUS**2

    # Number of particles needed to represent the dye (they overlap a little bit ~ 1.6 instead of 2)
    NUMBER_PARTICLES_BOTTOM_WALL = np.ceil(SL/(RADIUS*1.6))
    NUMBER_PARTICLES_SIDE_WALL = np.ceil(SH/(RADIUS*1.6))
    NUMBER_PARTICLES_DYE_BOTTOM_WALL = np.ceil(DL/(RADIUS*1.6))
    NUMBER_PARTICLES_DYE_SIDE_WALL = np.ceil(DH/(RADIUS*1.6))

    # Simulation parameters taken from BIRWISCH et al., 2009
    # Material
    N = 3000
    N += int(2*NUMBER_PARTICLES_BOTTOM_WALL + 2*NUMBER_PARTICLES_SIDE_WALL + NUMBER_PARTICLES_DYE_BOTTOM_WALL + 2*NUMBER_PARTICLES_DYE_SIDE_WALL) # Number of grains
    MASS= 7.63E-8 # Mass of each grain --- (kg)
    MU = 0.5 # Drag coefficient --- (dimensionless)
    MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
    KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
    MU_W = 0.15 # --- (dimensionless)

    # Repulsion parameters
    #E_TILDE = 1.0E-30 # Young's modulus / (1 - nu^2) --- (Pa)
    E_TILDE = 8.0E32 * 1.0E-30 # Young's modulus / (1 - nu^2) --- (Pa)
    #E_TILDE = 15.0E32 * 1.0E-30# For 5.E-4 radius
    #E_TILDE = 0.

    # Viscous parameter
    GAMMA_R = 2.0E-4 * 1.0E6 # Gamma / R --- (Pa . s/m)

    # Friction parameter
    #GBPM_GAMMA = 1.E-5
    GBPM_GAMMA = 4.E-5

    # Material of wall
    WALL_MASS = INFINITY

    # Misc parameters
    T0 = 0. # Initial time --- (s)
    #TF = 3. # End time --- (s)
    STEPS = 1000 # Number of steps --- (integer)
    #STEPS = 5 # Number of steps --- (integer)
    #DT = (TF-T0)/STEPS # Timestep between iterations --- (s)
    DT = 1.E-4/4.

    # Save parameters
    SAVE_ENABLED = False
    SAVE_SESSION_STEP_INTERVAL = 400 # Number of steps between saving session
    SAVE_SESSION_OUTPUT_PATH = "../output/simulation" + "_RADIUS" + str(RADIUS) + "_DT" + str(DT) + "_ETILDE" + str(E_TILDE) + "_GAMMAR" + str(GAMMA_R) + "_GBPMGAMMA" + str(GBPM_GAMMA) + "_N" + str(N)
    SAVE_SESSION_DIFFERENT_FILE_PER_STEP = True # Set if you want to progressively export the simulation state

    if SAVE_ENABLED and os.path.exists(SAVE_SESSION_OUTPUT_PATH):
        print("Output path already exists. Exiting...")
        exit()

# Derivative constants
TF = T0 + STEPS*DT
EFFECTIVE_RADIUS = (RADIUS*RADIUS)/(RADIUS+RADIUS)
