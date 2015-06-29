from parameters import *
from matrix_initialization import *
from forces import *

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
    plt.text(0.7*xf,0.87*yf, 'RADIUS = ' + str(RADIUS))
    plt.text(0.7*xf,0.84*yf, 'E_TILDE = ' + str(E_TILDE))
    plt.text(0.7*xf,0.81*yf, 'GAMMA_R = ' + str(GAMMA_R))
    plt.text(0.7*xf,0.78*yf, 'GBPM_GAMMA = ' + str(GBPM_GAMMA))
    plt.text(0.7*xf,0.75*yf, 'N = ' + str(N))

# Initial shoe velocity
shoe_velocity = 0.0

print 'Beggining solving steps...'

# Solving steps
for step in range(1, STEPS+1):

    time = T0 + (step-1)*DT

    # First, we sort by y position to optimize contact forces
    current_matrix = current_matrix[current_matrix[:,Y].argsort()]

    # Store current matrix as last matrix
    last_matrix = np.copy(current_matrix)

    # Reset matrix of forces
    current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

    current_matrix = apply_forces(current_matrix)

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

    print 'Done step ', str(step), '.' 

print 'Finished solving steps...'

