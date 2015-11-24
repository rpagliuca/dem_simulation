# -*- coding: utf-8 -*-

import os
import forces
import matplotlib.pyplot as plt
import parameters as p # Load all global variables from parameters
import numpy as np
import joblib

def main_loop(current_matrix):

    global shoe_velocity

    # Initializing matplotlib scatter plot data
    if p.realtimePlot or p.stepPlotFlag:

        scatterPointsParticles = None

        # Original limits
        print p.max_x
        print p.max_y
        x0 = p.min_x; xf = p.max_x; xl = xf-x0
        y0 = p.min_y; yf = p.max_y; yl = yf-y0

        # Limits with a little margin
        margin = 0.05
        x0 = x0 - margin*xl; xf = xf + margin*xl;
        y0 = y0 - margin*yl; yf = yf + margin*yl;

        # Initializate matplotlib pyplot
        plt.ion()
        plt.show()
        main_axes = plt.gca()
        main_axes.axis([x0, xf, y0, yf])
        main_axes.set_aspect('equal')

        # Initializate static plot text
        static_text = 'DT = ' + str(p.DT) + "\n"
        static_text += 'RADIUS = ' + str(p.RADIUS) + "\n"
        static_text += 'E_TILDE = ' + str(p.E_TILDE) + "\n"
        static_text += 'GAMMA_R = ' + str(p.GAMMA_R) + "\n"
        static_text += 'GBPM_GAMMA = ' + str(p.GBPM_GAMMA) + "\n"
        static_text += 'N_PARTICLES = ' + str(p.N_PARTICLES) + "\n"
        static_text += 'N_WALL = ' + str(p.N_WALL) + "\n"

    print 'Beggining solving steps...'

    # Converting current_matrix to memmap
    mem_share_folder = "/dev/shm"
    mem_share_file = os.path.join(mem_share_folder, "dem_simulation.mmap")
    # Create mem_share_folder if not exists
    if os.path.exists(mem_share_file):
        os.unlink(mem_share_file)
    joblib.dump(current_matrix, mem_share_file)
    current_matrix_ram = current_matrix
    current_matrix = joblib.load(mem_share_file, mmap_mode="r+")
    current_matrix[:,:] = current_matrix_ram[:,:]

    # Solving steps
    for step in range(p.start_step, p.STEPS+1):

        time = p.T0 + (step-1)*p.DT

        if p.simulation_mode != 'replay':

            # First, we sort by y position to optimize contact forces
            current_matrix[:,:] = current_matrix[current_matrix[:,p.Y].argsort()]

            # Store current matrix as last matrix
            last_matrix = np.copy(current_matrix)

            # Reset matrix of forces
            current_matrix[:,p.FX:p.FY+1] = np.zeros((p.N,p.DIMENSIONS))

            current_matrix[:,:] = forces.apply_forces(current_matrix)

            # Calculate velocities from force
            current_matrix[:,p.VX:p.VY+1] = last_matrix[:,p.VX:p.VY+1] + ((current_matrix[:,p.FX:p.FY+1] + last_matrix[:,p.FX:p.FY+1]).transpose()*current_matrix[:,p.T]/(2*current_matrix[:,p.M])).transpose() * p.DT

            # Move shoe horizontally
            if (time >= 0.2):
                p.shoe_velocity = 0.09
            current_matrix[current_matrix[:,p.WT] == 2, p.VX] = p.shoe_velocity

            # Calculate position from force and velocity
            current_matrix[:,p.X:p.Y+1] = last_matrix[:,p.X:p.Y+1] + current_matrix[:,p.VX:p.VY+1]*p.DT + (current_matrix[:,p.FX:p.FY+1].transpose()/(2*current_matrix[:,p.M])).transpose() * p.DT**2

        if p.realtimePlot or (p.stepPlotFlag and step % p.stepPlotSteps == 0):
            plotNow = True
        else:
            plotNow = False

        if plotNow:
            if scatterPointsParticles:
                scatterPointsParticles.remove()
                scatterPointsWall1.remove()
                scatterPointsWall2.remove()
                text.remove()
            scatterPointsParticles = plt.scatter(current_matrix[current_matrix[:,p.T] == 1, p.X], current_matrix[current_matrix[:, p.T] == 1, p.Y], s=p.scatterPlotPointSize, facecolors='none')
            scatterPointsWall1 = plt.scatter(current_matrix[current_matrix[:, p.WT] == 1, p.X], current_matrix[current_matrix[:, p.WT] == 1, p.Y], s=p.scatterPlotPointSize, facecolors='none', color='green')
            scatterPointsWall2 = plt.scatter(current_matrix[current_matrix[:, p.WT] == 2, p.X], current_matrix[current_matrix[:, p.WT] == 2, p.Y], s=p.scatterPlotPointSize, facecolors='none', color='red')

            # Custom function to export replay plots as pdfs
            SAVE_PLOT = False
            if SAVE_PLOT:
                plt.gcf().set_size_inches(10, 7)
                text = main_axes.text(0.95, 0.05, '#' + str(step), horizontalalignment='right', verticalalignment='bottom', size='20', transform=main_axes.transAxes)
                plt.xticks([])
                plt.yticks([])
                output_path = p.SAVE_SESSION_OUTPUT_PATH
                output_image_dir = os.path.join(output_path, "..", "..", "replay_images")
                output_image = os.path.join(output_image_dir, 'step' + str(step).zfill(11) + '.pdf')
                print "Saving image " + output_image + "..."
                if not os.path.exists(output_image_dir):
                    os.makedirs(output_image_dir)
                plt.savefig(output_image, format='pdf', bbox_inches='tight', dpi=50)

            else:
                time_text = 'Time: ' + format(time, '.5f') + 's of ' + format(p.TF, '.5f') + "s\n"
                text = main_axes.text(0.95, 0.95, time_text + static_text, horizontalalignment='right', verticalalignment='top', transform=main_axes.transAxes)
                plt.draw()
                plt.pause(1.0E-1)

        if p.SAVE_ENABLED and step % p.SAVE_SESSION_STEP_INTERVAL == 0:
            output_path = p.SAVE_SESSION_OUTPUT_PATH
            if (p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP):
                output_path = os.path.join(output_path, "/step" + str(step))

            print "Saving session to folder " + output_path + "/ ..."
            # Create output folder if not exists
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # Unify all parameters into a single numpy array
            parameters = np.array([p.SAVE_SESSION_STEP_INTERVAL, p.SAVE_SESSION_DIFFERENT_FILE_PER_STEP, p.SH, p.SL, p.SH_MULTIPLICATOR, p.DH, p.DL, p.N, p.RADIUS, p.scatterPlotPointSize, p.MASS, p.MU, p.MU_A, p.KAPPA_R, p.MU_W, p.shoe_velocity, step, p.T0, p.DT, p.STEPS, p.GBPM_GAMMA, p.GAMMA_R, p.E_TILDE, p.N_WALL, p.N_PARTICLES, p.min_x, p.min_y, p.max_x, p.max_y ])

            # Export current state to file
            np.savetxt(os.path.join(output_path, "/current_matrix.txt"), current_matrix)
            # Export simulation parameters to file
            np.savetxt(os.path.join(output_path, "/parameters.txt"), parameters)

        print 'Done step ', str(step), '.' 

    plt.clf()
    
    print 'Finished solving steps...'
