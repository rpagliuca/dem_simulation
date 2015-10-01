#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os
import parameters as p
import glob
import main_loop
import saved_states
import dye_profiles.h as dye_profile
import matplotlib.pyplot as plt

def post_analysis_density(import_path):
    # 1) Load state
    saved_states.load_state(import_path)

    # 2) Create mesh having cells of 0.1 diam x 0.1 diam in the region where we have powder particles (!excluding wall particles!)
    x_min = p.current_matrix[p.current_matrix[:, p.T] == 1, p.X].min()
    x_max = p.current_matrix[p.current_matrix[:, p.T] == 1, p.X].max()
    y_min = p.current_matrix[p.current_matrix[:, p.T] == 1, p.Y].min()
    y_max = p.current_matrix[p.current_matrix[:, p.T] == 1, p.Y].max()

    # Get actual region of the dye, to ignore distant particles
    dye_profile.init()
    dye_profile.draw([], True)
    dye_region = dye_profile.get_dye_region()

    # Ignore points if outside dye area
    #x_min = max(x_min, dye_region['x_min'])
    #x_max = min(x_max, dye_region['x_max'])
    #y_min = max(y_min, dye_region['y_min'])

    #x_min = 0.051
    #x_max = 0.053
    #y_min = -0.014

    cell_width = cell_height = 0.5 * p.RADIUS*2.0
    num_x = int(np.ceil((x_max - x_min)/cell_width))
    num_y = int(np.ceil((y_max - y_min)/cell_height))

    # 3) Loop over every mesh cell
        # 3.1 ) Loop over every particle and sum the mesh cell area for the nearest particle

    # First we zero the count of owned cells

    # I'm using columns VY and VX, not needed otherwise in this function, to count the total of mesh cells "owned" by a particle and the distance to the current cell
    p.current_matrix[:,p.VY] = 0
    p.current_matrix[p.current_matrix[:, p.T] != 1, p.VX] = p.INFINITY
    for i in range (0, num_x):
        for j in range(0, num_y):
            cell_x = x_min + i*cell_width
            cell_y = y_min + j*cell_height
            distances = np.square(p.current_matrix[:, p.X:p.Y+1] - np.matrix([cell_x, cell_y]))
            p.current_matrix[:, p.VX] = np.sqrt(distances[:, 0] + distances[:, 1]).transpose()
            # Now we sort by distance (here called VX)
            p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VX].argsort()]
            # The first particle after sorting will have an increment on the number of owned cells
            p.current_matrix[0,p.VY] = p.current_matrix[0,p.VY] + 1

        # Show progress of calculation
        print float(round(float(i)/float(num_x)*1000))/10.0, "%"

    # 4) Plot intensity
    p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VY].argsort()]
    plt.show()
    # Limits with a little margin
    margin = 0.05
    xl = x_max - x_min; yl = y_max - y_min
    x0 = x_min - margin*xl; xf = x_max + margin*xl;
    y0 = y_min - margin*yl; yf = y_max + margin*yl;
    plt.axis([x0, xf, y0, yf])
    plt.axes().set_aspect('equal')
    sc = plt.scatter(p.current_matrix[p.current_matrix[:,p.T] == 1, p.X], p.current_matrix[p.current_matrix[:, p.T] == 1, p.Y], c=np.pi*(p.RADIUS)**2.0/(p.current_matrix[p.current_matrix[:, p.T] == 1, p.VY]*cell_width*cell_height), s=50, linewidth=0, cmap=plt.cm.summer_r)
    plt.colorbar(sc)
    plt.show()
