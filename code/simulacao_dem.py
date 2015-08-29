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
    saved_states.load_state(p.saved_state_path)

# Calculate some derivatives parameters
p.load_parameters_post()

# Run the simulation per se
main_loop.main_loop(p.current_matrix)
