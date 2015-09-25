#!/usr/bin/python
#-*- coding: utf-8 -*-

import numpy as np
import os
import parameters as p
import glob
import main_loop
import saved_states

def post_analysis_density(import_path):
    # 1) Load state
    saved_states.load_state(import_path)

    # 2) Create mesh having cells of 0.1 diam x 0.1 diam in the region where we have powder particles (!excluding wall particles!)
    x_min = p.current_matrix[p.current_matrix[:, p.T] == 1, p.X].min()
    x_max = p.current_matrix[p.current_matrix[:, p.T] == 1, p.X].max()
    y_min = p.current_matrix[p.current_matrix[:, p.T] == 1, p.Y].min()
    y_max = p.current_matrix[p.current_matrix[:, p.T] == 1, p.Y].max()

    cell_width = cell_height = 0.1 * p.RADIUS*2.0
    num_x = int(np.ceil((x_max - x_min)/cell_width))
    num_y = int(np.ceil((y_max - y_min)/cell_height))
    
    # 3) Loop over every mesh cell
        # 3.1 ) Loop over every particle and sum the mesh cell area for the nearest particle

    # First we zero the count of owned cells
    p.current_matrix[:,p.VY] = 0
    p.current_matrix[p.current_matrix[:, p.T] != 1, p.VX] = p.INFINITY
    for i in range (0, num_x):
        for j in range(0, num_y):
            cell_x = x_min + i*cell_width
            cell_y = y_min + j*cell_height
            # I'm using column VX and VY, irrelevant here, to count the total of mesh cells "owned" by a particle and the current distance
            distances = np.square(p.current_matrix[p.current_matrix[:, p.T] == 1, p.X:p.Y+1] - np.array(cell_x, cell_y))
            p.current_matrix[p.current_matrix[:, p.T] == 1, p.VX] = np.sqrt(distances[:, 0] + distances[:, 1])
            # Now we sort by distance (here called VX)
            p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VX].argsort()]
            # The first particle after sorting will have an increment on the number of owned cells
            p.current_matrix[0,p.VY] += 1
            print p.current_matrix[0, p.VY]

        # Show progress of calculation
        print float(round(float(i)/float(num_x)*1000))/10.0, "%"
