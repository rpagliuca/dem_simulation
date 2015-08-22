#!/usr/bin/python
# -*- coding: utf-8 -*-

import matrix_initialization
import main_loop
from parameters import * # Load all global variables from parameters

if not load_simulation:
    current_matrix = matrix_initialization.matrix_initialization()

main_loop.main_loop(current_matrix)
