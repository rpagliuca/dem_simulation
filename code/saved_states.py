#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os

def saved_states():

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
