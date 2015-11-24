#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import matrix_initialization
import time
import os
import sys
import datetime

def load_parameters_pre():

    global min_x, min_y, max_x, max_y, DESIRED_N_PARTICLES, simulation_mode, stepPlotFlag, stepPlotSteps, number_of_cores, saved_state_path, realtimePlot, stepPlotFlag, lastPlotFlag, G, PI, INFINITY, X, Y, VX, VY, FX, FY, M, T, WT, DIMENSIONS

    # Reproduce random results for debugging
    np.random.seed(1)

    # Multicore
    number_of_cores = 8 # Used if multicore/multithreaded simulation enabled on forces.py

    # Read some parameters from command line
    if len(sys.argv) >= 3:
        # path of the saved state, Used for both 'load' and 'replay' modes
        saved_state_path = sys.argv[2]
    if len(sys.argv) >= 2:
        # simulation_mode can be 'new', 'load' or 'replay'
        simulation_mode = sys.argv[1]
    else:
        print "Usage: python " + sys.argv[0] + " [new|load|replay|list] <import_path>"
        exit()

    # Default flags
    realtimePlot = False # realtimePlot will slow down simulation, but gives a nice feedback
    stepPlotFlag = True # Use with stepPlotSteps to plot every X steps
    stepPlotSteps = 20

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

    min_x = 0
    max_x = 0
    min_y = 0
    max_y = 0

# Default parameters for new simulations
def load_default_parameters():

    global DESIRED_N_PARTICLES, shoe_velocity, start_step, SH, SH_MULTIPLICATOR, SL, DH, DL, RADIUS, scatterPlotPointSize, NUMBER_PARTICLES_BOTTOM_WALL, NUMBER_PARTICLES_SIDE_WALL, NUMBER_PARTICLES_DYE_BOTTOM_WALL, NUMBER_PARTICLES_DYE_SIDE_WALL, N, MASS, MU, MU_A, KAPPA_R, MU_W, E_TILDE, GAMMA_R, GBPM_GAMMA, WALL_MASS, T0, STEPS, DT, SAVE_ENABLED, SAVE_SESSION_STEP_INTERVAL, SAVE_SESSION_OUTPUT_PATH, SAVE_SESSION_DIFFERENT_FILE_PER_STEP

    # Initial shoe velocity
    shoe_velocity = 0.0

    # Start on step 1 for new simulations
    start_step = 1

    # Particle size
    RADIUS = 1.5E-4 # Radius of each grain --- (m)
    scatterPlotPointSize = 1.0E8 * RADIUS**2

    # This is the desired number of particles to be simulated
    DESIRED_N_PARTICLES = 10000
    # Note: the real number of simulated particles (N) will vary. It is the sum of the DESIRED_N_PARTICLES, with the number of particles needed to create the walls of the dye, and the subtraction of the conflicting overlapped particles that will be excluded from simulation after random generation

    # Simulation parameters taken from BIRWISCH et al., 2009

    # Material
    MASS= 7.63E-8 # Mass of each grain --- (kg)
    MU = 0.5 # Drag coefficient --- (dimensionless)
    MU_A = 2.0E-5 # Stoke air drag coefficient --- (Pa . s)
    KAPPA_R = 1.0E6 # Kappa / R --- (Pa)
    MU_W = 0.15 # --- (dimensionless)

    # Repulsion parameters
    E_TILDE = 1.0E7 # Young's modulus / (1 - nu^2) --- (Pa)

    # Viscous parameter
    GAMMA_R = 1.0E6 # Gamma / R --- (Pa . s/m)

    # Friction parameter - this parameters is the only one which is not available on Bierwisch et al (2009),
    # because I have not modelled the friction force force as Cundall, but as Haff and Werner
    GBPM_GAMMA = 4.E-5

    # Material of wall
    WALL_MASS = INFINITY

    # Misc parameters
    T0 = 0. # Initial time --- (s)
    STEPS = 2000000 # Number of steps --- (integer)
    DT = 1.9E-6 

    N = 0 # Updated on matrix_initialization.py to hold the total number of particles

    # Save parameters
    SAVE_ENABLED = True
    SAVE_SESSION_STEP_INTERVAL = 1000 # Number of steps between saving session
    SAVE_SESSION_OUTPUT_PATH = "../output/"
    SAVE_SESSION_DIFFERENT_FILE_PER_STEP = True # Set if you want to progressively export the simulation state

def load_parameters_post():

    global simulation_mode, TF, EFFECTIVE_RADIUS, SAVE_SESSION_OUTPUT_PATH

    # Derivative constants
    TF = T0 + STEPS*DT
    EFFECTIVE_RADIUS = (RADIUS*RADIUS)/(RADIUS+RADIUS)

    # Define output filename of simulation
    if simulation_mode == 'new':
        data = datetime.datetime.now().isoformat()
        data = data[0:19].replace(':', '').replace('-', '')
        filename =  data + "_simulation_RADIUS" + str(RADIUS) + "_DT" + str(DT) + "_ETILDE" + str(E_TILDE) + "_GAMMAR" + str(GAMMA_R) + "_GBPMGAMMA" + str(GBPM_GAMMA) + "_N" + str(N)
        SAVE_SESSION_OUTPUT_PATH = os.path.join(SAVE_SESSION_OUTPUT_PATH, filename)
        if SAVE_ENABLED:
            if os.path.exists(SAVE_SESSION_OUTPUT_PATH):
                print("Output path already exists. Exiting...")
                exit()
