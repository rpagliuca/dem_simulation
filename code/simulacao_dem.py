# -*- coding: utf-8 -*-

import matrix_initialization
import main_loop
import saved_states
import parameters as p

# Function to list saved states
#saved_states.saved_states()

# Define some basic constants and parameters shared between new and saved simulations
p.load_parameters_pre()

# Load parameters differently, depending if it was a new simulation, or loading from a saved state
if not p.load_saved_state:
    # Load default parameters
    p.load_default_parameters()
    matrix_initialization.matrix_initialization()
else: # Or load parameters from saved state
    saved_states.load_state("/home/rsantos/Desktop/simulacao_dem/output/simulation_RADIUS0.0008_DT0.0001_ETILDE800.0_GAMMAR200.0_GBPMGAMMA4e-05_N3208/step800")

# Calculate some derivatives parameters
p.load_parameters_post()

main_loop.main_loop(p.current_matrix)
