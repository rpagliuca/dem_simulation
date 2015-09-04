#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os
import parameters as p
import glob
import main_loop

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
        p.N = int(p.N)

        p.start_step = p.step + 1


def replay(import_path):

    print "Replaying simulation from:"
    print import_path
    print ""

    # Search and return step* folder inside import_path (e.g. step100, step200, step300 etc), ordered by number
    steps_paths = sorted(glob.glob(os.path.join(import_path, "step*")), key=lambda name: int(name.replace(import_path, "").replace("step", "")))

    if len(steps_paths) == 0:
        print "There is no step* folder to replay from."
        exit()
    else:
        for step_path in steps_paths:
            load_state(os.path.join(step_path))
            p.load_parameters_post()
            p.STEPS = p.start_step
            p.realtimePlot = True
            main_loop.main_loop(p.current_matrix)

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
