from parameters import *
from matrix_initialization import *
from functions import *

# Initializing matplotlib scatter plot data
stepPlot = False
lastPlot = False
if realtimePlot:
    scatterPoints = None
    plt.ion()
    plt.show()
    plt.axis([0 - L*0.05, L + L*0.05, 0 - 2*H*0.05, 2*H + 2*H*0.05])
    plt.axes().set_aspect('equal')

# Solving steps
for step in range(1, STEPS+1):

    # First, we sort by y position to optimize contact forces
    current_matrix = current_matrix[current_matrix[:,Y].argsort()]

    # Store current matrix as last matrix
    last_matrix = np.copy(current_matrix)

    # Reset matrix of forces
    current_matrix[:,FX:FY+1] = np.zeros((N,DIMENSIONS))

    current_matrix = apply_forces(current_matrix)

    # Calculate velocities from force
    current_matrix[:,VX:VY+1] = last_matrix[:,VX:VY+1] + ((current_matrix[:,FX:FY+1] + last_matrix[:,FX:FY+1]).transpose()*current_matrix[:,T]/(2*current_matrix[:,M])).transpose() * DT

    # Calculate position from force and velocity
    #current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + last_matrix[:,VX:VY+1]*DT + (last_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2
    current_matrix[:,X:Y+1] = last_matrix[:,X:Y+1] + current_matrix[:,VX:VY+1]*DT + (current_matrix[:,FX:FY+1].transpose()/(2*current_matrix[:,M])).transpose() * DT**2

    # Draw animation frame
    if step == STEPS:
        lastPlot = True
    if stepPlotFlag and (step == STEPS or step == np.round(STEPS/2) or step == np.round(STEPS/4)):
        stepPlot = True
    if not realtimePlot and (stepPlot or (lastPlotFlag and lastPlot)):
        scatterPoints = None
        plt.ion()
        plt.show()
        plt.axis([0 - L*0.05, L + L*0.05, 0 - 2*H*0.05, 2*H + 2*H*0.05])
        plt.axes().set_aspect('equal')
        realtimePlot = True

    if realtimePlot:
        if scatterPoints:
            scatterPoints.remove()
        scatterPoints = plt.scatter(current_matrix[:, X], current_matrix[:,Y], s=scatterPlotPointSize, facecolors='none')
        if stepPlot or lastPlot:
            plt.ioff()
            plt.show()
        else:
            plt.draw()
            plt.pause(1.0E-10)

    if stepPlot:
        realtimePlot = False
        stepPlot = False

    print 'Done step ', str(step), '.' 

