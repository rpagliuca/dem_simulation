# -*- coding: utf-8 -*-

import os
import forces
import matplotlib.pyplot as plt
import parameters as p # Load all global variables from parameters
import numpy as np

def main_loop(current_matrix):

    global shoe_velocity

# Initializing matplotlib scatter plot data
    if p.realtimePlot or p.stepPlotFlag:
        scatterPoints = None
        plt.ion()
        plt.show()
        # Original limits
        x0 = 0; xf = 2*p.SL+p.DL; xl = xf-x0
        y0 = -p.DH; yf = p.SH_MULTIPLICATOR * p.SH; yl = yf-y0
        # Limits with a little margin
        margin = 0.05
        x0 = x0 - margin*xl; xf = xf + margin*xl;
        y0 = y0 - margin*yl; yf = yf + margin*yl;
        plt.axis([x0, xf, y0, yf])
        plt.axes().set_aspect('equal')
        plt.text(0.7*xf,0.87*yf, 'DT = ' + str(p.DT))
        plt.text(0.7*xf,0.84*yf, 'RADIUS = ' + str(p.RADIUS))
        plt.text(0.7*xf,0.81*yf, 'E_TILDE = ' + str(p.E_TILDE))
        plt.text(0.7*xf,0.78*yf, 'GAMMA_R = ' + str(p.GAMMA_R))
        plt.text(0.7*xf,0.75*yf, 'GBPM_GAMMA = ' + str(p.GBPM_GAMMA))
        plt.text(0.7*xf,0.72*yf, 'N = ' + str(p.N))

    print 'Beggining solving steps...'

    # Solving steps
    for step in range(p.start_step, p.STEPS+1):

        time = p.T0 + (step-1)*p.DT

        # First, we sort by y position to optimize contact forces
        current_matrix = current_matrix[current_matrix[:,p.Y].argsort()]

        # Store current matrix as last matrix
        last_matrix = np.copy(current_matrix)

        # Reset matrix of forces
        current_matrix[:,p.FX:p.FY+1] = np.zeros((p.N,p.DIMENSIONS))

        #current_matrix = forces.apply_forces(current_matrix)
        forces.apply_forces(current_matrix)

        # Calculate velocities from force
        current_matrix[:,p.VX:p.VY+1] = last_matrix[:,p.VX:p.VY+1] + ((current_matrix[:,p.FX:p.FY+1] + last_matrix[:,p.FX:p.FY+1]).transpose()*current_matrix[:,p.T]/(2*current_matrix[:,p.M])).transpose() * p.DT

        # Move shoe horizontally
        if (time == 0.1):
            p.shoe_velocity = 0.2
        current_matrix[current_matrix[:,p.WT] == 2, p.VX] = p.shoe_velocity

        # Calculate position from force and velocity
        #current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + (last_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2
        current_matrix[:,p.X:p.Y+1] = last_matrix[:,p.X:p.Y+1] + current_matrix[:,p.VX:p.VY+1]*p.DT + (current_matrix[:,p.FX:p.FY+1].transpose()/(2*current_matrix[:,p.M])).transpose() * p.DT**2

        if p.realtimePlot or (p.stepPlotFlag and step % p.stepPlotSteps == 0):
            plotNow = True
        else:
            plotNow = False

        if plotNow:
            if scatterPoints:
                scatterPoints.remove()
                text.remove()
            scatterPoints = plt.scatter(current_matrix[:, p.X], current_matrix[:,p.Y], s=p.scatterPlotPointSize, facecolors='none')
            text = plt.text(0.7*xf,0.9*yf, 'Time: ' + format(time, '.5f') + 's of ' + format(p.TF, '.5f') + 's')
            plt.draw()
            plt.pause(1.0E-10)

        if p.SAVE_ENABLED and step % p.SAVE_SESSION_STEP_INTERVAL == 0:
            output_path = p.SAVE_SESSION_OUTPUT_PATH
            if (p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP):
                output_path = output_path + "/step" + str(step)

            print "Saving session to folder " + output_path + "/ ..."
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # Unify all parameters into a numpy array
            parameters = np.array([p.SAVE_SESSION_STEP_INTERVAL, p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP, p.SH, p.SL, p.SH_MULTIPLICATOR, p.DH, p.DL, p.N, p.RADIUS, p.scatterPlotPointSize, p.MASS, p.MU, p.MU_A, p.KAPPA_R, p.MU_W, p.shoe_velocity, step, p.T0, p.DT, p.STEPS, p.GBPM_GAMMA, p.GAMMA_R, p.E_TILDE])

            np.savetxt(output_path + "/current_matrix.txt", current_matrix)
            np.savetxt(output_path + "/parameters.txt", parameters)

        print 'Done step ', str(step), '.' 

    print 'Finished solving steps...'
