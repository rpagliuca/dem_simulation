#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import time

from parameters import *
from matrix_initialization import *
from forces import *
from main_loop import *

if not load_simulation:
    current_matrix = matrix_initialization()

main_loop(current_matrix)
