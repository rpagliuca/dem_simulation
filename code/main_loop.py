import os
import forces
import matplotlib.pyplot as plt
from parameters import * # Load all global variables from parameters

def main_loop(current_matrix):

    global shoe_velocity

# Initializing matplotlib scatter plot data
    stepPlot = False
    if realtimePlot:
        scatterPoints = None
        plt.ion()
        plt.show()
        # Original limits
        x0 = 0; xf = 2*SL+DL; xl = xf-x0
        y0 = -DH; yf = SH_MULTIPLICATOR * SH; yl = yf-y0
        # Limits with a little margin
        margin = 0.05
        x0 = x0 - margin*xl; xf = xf + margin*xl;
        y0 = y0 - margin*yl; yf = yf + margin*yl;
        plt.axis([x0, xf, y0, yf])
        plt.axes().set_aspect('equal')
        plt.text(0.7*xf,0.87*yf, 'DT = ' + str(DT))
        plt.text(0.7*xf,0.84*yf, 'RADIUS = ' + str(RADIUS))
        plt.text(0.7*xf,0.81*yf, 'E_TILDE = ' + str(E_TILDE))
        plt.text(0.7*xf,0.78*yf, 'GAMMA_R = ' + str(GAMMA_R))
        plt.text(0.7*xf,0.75*yf, 'GBPM_GAMMA = ' + str(GBPM_GAMMA))
        plt.text(0.7*xf,0.72*yf, 'N = ' + str(N))

    print 'Beggining solving steps...'

    # Solving steps
    for step in range(start_step, STEPS+1):

        time = T0 + (step-1)*DT

        # First, we sort by y position to optimize contact forces
        current_matrix = current_matrix[current_matrix[:,Y].argsort()]

        # Store current matrix as last matrix
        last_matrix = np.copy(current_matrix)

        # Reset matrix of forces
        current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

        #current_matrix = forces.apply_forces(current_matrix)
        forces.apply_forces(current_matrix)

        # Calculate velocities from force
        current_matrix[:,VX:VY+1] = last_matrix[:,VX:VY+1] + ((current_matrix[:,FX:FY+1] + last_matrix[:,FX:FY+1]).transpose()*current_matrix[:,T]/(2*current_matrix[:,M])).transpose() * DT

        # Move shoe horizontally
        if (time == 0.1):
            shoe_velocity = 0.2
        current_matrix[current_matrix[:,WT] == 2, VX] = shoe_velocity

        # Calculate position from force and velocity
        #current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + (last_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2
        current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + current_matrix[:,VX:VY+1]*DT + (current_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2

        if realtimePlot:
            if scatterPoints:
                scatterPoints.remove()
                text.remove()
            scatterPoints = plt.scatter(current_matrix[:, X], current_matrix[:,Y], s=scatterPlotPointSize, facecolors='none')
            text = plt.text(0.7*xf,0.9*yf, 'Time: ' + format(time, '.5f') + 's of ' + format(TF, '.5f') + 's')
            plt.draw()
            plt.pause(1.0E-10)

        if SAVE_ENABLED and step % SAVE_SESSION_STEP_INTERVAL == 0:
            output_path = SAVE_SESSION_OUTPUT_PATH
            if (SAVE_SESSION_DIFFERENT_FILE_PER_STEP):
                output_path = output_path + "/step" + str(step)

            print "Saving session to folder " + output_path + "/ ..."
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # Unify all parameters into a numpy array
            parameters = np.array([SAVE_SESSION_STEP_INTERVAL, SAVE_SESSION_DIFFERENT_FILE_PER_STEP, SH, SL, SH_MULTIPLICATOR, DH, DL, N, RADIUS, scatterPlotPointSize, MASS, MU, MU_A, KAPPA_R, MU_W, shoe_velocity, step, T0, DT, STEPS, GBPM_GAMMA, GAMMA_R, E_TILDE])

            np.savetxt(output_path + "/current_matrix.txt", current_matrix)
            np.savetxt(output_path + "/parameters.txt", parameters)

        print 'Done step ', str(step), '.' 

    print 'Finished solving steps...'
