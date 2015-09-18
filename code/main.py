# -*- coding: utf-8 -*-

import matrix_initialization
import main_loop
import saved_states
import parameters as p

# Define some basic constants and parameters shared between new and saved simulations
p.load_parameters_pre()

# Function to list saved states

if p.simulation_mode == 'list':

    saved_states.list_saved_states(p.saved_state_path)

elif p.simulation_mode == 'new':

    # Load default parameters
    p.load_default_parameters()
    matrix_initialization.matrix_initialization()
    # Calculate some derivatives parameters
    p.load_parameters_post()
    # Run the simulation per se
    main_loop.main_loop(p.current_matrix)

elif p.simulation_mode == 'load': # Or load parameters from saved state

    saved_states.load_state(p.saved_state_path)
    # Calculate some derivatives parameters
    p.load_parameters_post()
    main_loop.main_loop(p.current_matrix)

elif p.simulation_mode == 'replay':

    # Replay simulation
    saved_states.replay(p.saved_state_path)
