# -*- coding: utf-8 -*-

import matrix_initialization
import main_loop
import saved_states
import parameters as p

#saved_states.saved_states()

if not load_simulation:
    current_matrix = matrix_initialization.matrix_initialization()

main_loop.main_loop(current_matrix)
