#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os
import parameters as p

def load_state(import_path):

    print "Loading state from:"
    print import_path
    print ""

    current_matrix_path = os.path.join(import_path, "current_matrix.txt")
    parameters_path = os.path.join(import_path, "parameters.txt")
    if not os.path.isfile(current_matrix_path) or not os.path.isfile(parameters_path):
        print "Error loading state from " + import_path + "."
        exit()
    else:
        p.current_matrix = np.loadtxt(current_matrix_path)
        p.SAVE_SESSION_STEP_INTERVAL, p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP, p.SH, p.SL, p.SH_MULTIPLICATOR, p.DH, p.DL, p.N, p.RADIUS, p.scatterPlotPointSize, p.MASS, p.MU, p.MU_A, p.KAPPA_R, p.MU_W, p.shoe_velocity, p.step, p.T0, p.DT, p.STEPS, p.GBPM_GAMMA, p.GAMMA_R, p.E_TILDE = np.loadtxt(parameters_path, unpack = True)

        print "N = " + str(p.N)
        p.SAVE_ENABLED = False
        # Convert some floats to int

        p.step = int(p.step)
        p.STEPS = int(p.STEPS)
        p.STEPS = 1000
        p.N = int(p.N)

        p.start_step = p.step + 1

def list_saved_states():

    saved_states_path = "/home/rsantos/Desktop/simulacao_dem/output/"
    print "Listing saved states on " + saved_states_path + ":"

    for root, directories, filenames in os.walk(saved_states_path):
        for directory in directories:
            import_path = os.path.join(root, directory) 
            clean_import_path = import_path.replace(saved_states_path, "")
            parameters_path = os.path.join(import_path, "parameters.txt")
            if os.path.isfile(parameters_path):
                SAVE_SESSION_STEP_INTERVAL, SAVE_SESSION_DIFFERENT_FILE_PER_STEP, SH, SL, SH_MULTIPLICATOR, DH, DL, N, RADIUS, scatterPlotPointSize, MASS, MU, MU_A, KAPPA_R, MU_W, shoe_velocity, step, T0, DT, STEPS, GBPM_GAMMA, GAMMA_R, E_TILDE = np.loadtxt(parameters_path, unpack = True)

                print ""
                print "==================================="
                print "State " + clean_import_path
                print "==================================="
                print "N = " + str(N)
                print "step = " + str(step)
                print "RADIUS = " + str(RADIUS)
                print "T0 = " + str(T0)
                print "DT = " + str(DT)
                print "GBPM_GAMMA = " + str(GBPM_GAMMA)
                print "GAMMA_R = " + str(GAMMA_R)
                print "E_TILDE = " + str(E_TILDE)
                print "MASS = " + str(MASS)
                print "MU = " + str(MU)
                print "MU_A = " + str(MU_A)
                print "KAPPA_R = " + str(KAPPA_R)
