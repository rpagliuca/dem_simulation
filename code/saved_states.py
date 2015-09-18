#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os
import parameters as p
import glob
import main_loop

def load_state(import_path, child_import = False, print_messages = True):

    if print_messages:
        print "Loading state from:"
        print import_path
        print ""

    if not child_import:
        p.SAVE_SESSION_OUTPUT_PATH = import_path

    current_matrix_path = os.path.join(import_path, "current_matrix.txt")
    parameters_path = os.path.join(import_path, "parameters.txt")

    if not os.path.isfile(current_matrix_path) or not os.path.isfile(parameters_path):
        steps_paths = glob.glob(os.path.join(import_path, "step*"))
        if len(steps_paths) == 0:
            if print_messages:
                print "Error loading state. Files current_matrix.txt and parameters.txt not found."
            return False
        else:
            steps_paths = sorted(steps_paths, key=lambda name: int(name.replace(import_path, "").replace("step", "").replace("/", "")))
            return load_state(steps_paths[-1], child_import = True, print_messages = print_messages)
    else:
        p.current_matrix = np.loadtxt(current_matrix_path)
        try:
            p.SAVE_SESSION_STEP_INTERVAL, p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP, p.SH, p.SL, p.SH_MULTIPLICATOR, p.DH, p.DL, p.N, p.RADIUS, p.scatterPlotPointSize, p.MASS, p.MU, p.MU_A, p.KAPPA_R, p.MU_W, p.shoe_velocity, p.step, p.T0, p.DT, p.STEPS, p.GBPM_GAMMA, p.GAMMA_R, p.E_TILDE, p.N_WALL, p.N_PARTICLES = np.loadtxt(parameters_path, unpack = True)
        except:
            if print_messages:
                print "Error: Unknown format for parameter.txt file."
            return False
        p.SAVE_ENABLED = True
        # Convert some floats to int
        p.step = int(p.step)
        p.STEPS = 1000000
        p.N = int(p.N)
        p.N_PARTICLES = p.N
        p.N_WALL = p.N
        p.start_step = p.step + 1
        return True

def replay(import_path):

    print "Replaying simulation from:"
    print import_path
    print ""

    # Search and return step* folder inside import_path (e.g. step100, step200, step300 etc), ordered by number
    steps_paths = glob.glob(os.path.join(import_path, "step*"))

    if len(steps_paths) == 0:
        print "There is no step* folder to replay from."
        exit()
    else:
        steps_paths = sorted(steps_paths, key=lambda name: int(name.replace(import_path, "").replace("step", "").replace("/", "")))
        step_number = 0
        for step_path in steps_paths:
            step_number += 1
            if step_number % 10 == 0: # Load every X step
                load_state(os.path.join(step_path))
                p.load_parameters_post()
                p.STEPS = p.start_step
                p.realtimePlot = True
                main_loop.main_loop(p.current_matrix)

def list_saved_states(saved_states_path):

    print "Listing saved states on " + saved_states_path + ":"

    for directory in glob.glob(os.path.join(saved_states_path, '*')):
        if load_state(directory, print_messages = False):
            print ""
            print "==================================="
            print "State " + directory
            print "==================================="
            print "N = " + str(p.N)
            print "step = " + str(p.step)
            print "RADIUS = " + str(p.RADIUS)
            print "T0 = " + str(p.T0)
            print "DT = " + str(p.DT)
            print "GBPM_GAMMA = " + str(p.GBPM_GAMMA)
            print "GAMMA_R = " + str(p.GAMMA_R)
            print "E_TILDE = " + str(p.E_TILDE)
            print "MASS = " + str(p.MASS)
            print "MU = " + str(p.MU)
            print "MU_A = " + str(p.MU_A)
            print "KAPPA_R = " + str(p.KAPPA_R)
        else:
            print ""
            print "==================================="
            print "Error loading " + directory
            print "==================================="
