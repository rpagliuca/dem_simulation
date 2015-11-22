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
import time
import threading
import datetime

density_polygon_list = {'x': [], 'y': []} # X, Y
density_scatter_points = None
density_scatter_points_lines = None

def plot_onclick(event):
    global density_polygon_list, density_scatter_points
    if event.dblclick:
        density_polygon_list['x'].append(event.xdata)
        density_polygon_list['y'].append(event.ydata)
        density_scatter_points_lines = plt.plot(density_polygon_list['x'], density_polygon_list['y'], '-o', c = 'lime')
        plt.draw()

def plot_keyrelease(event):
    global density_polygon_list, density_scatter_points
    if event.key == 'enter':
        density_polygon_list['x'].append(density_polygon_list['x'][0])
        density_polygon_list['y'].append(density_polygon_list['y'][0])
        density_scatter_points_lines = plt.plot(density_polygon_list['x'], density_polygon_list['y'], '-o', c = 'lime')
        #plt.draw()
        plt.close()
        # I commented the threading, because it was throwing random errors
        #timer = threading.Timer(2, lambda: plt.close()).start()

def point_belongs_to_polygon(x0, y0, polygon_ordered_list):
    p = polygon_ordered_list
    # Generate list of y = a*x + b for every polygon pair 
    point_count_on_the_right = 0
    for i in range(0, len(polygon_ordered_list['x'])-1):
        try:
            a = (polygon_ordered_list['y'][i+1] - polygon_ordered_list['y'][i]) / (polygon_ordered_list['x'][i+1] - polygon_ordered_list['x'][i])
        except:
            # Division by zero, it means the points have the same x (vertical line)
            if polygon_ordered_list['x'][i] > x0:
                point_count_on_the_right += 1
            continue

        # If not division by zero, continue
        b = polygon_ordered_list['y'][i] - a * polygon_ordered_list['x'][i]
        x = lambda y: (y-b)/a
        if min(p['y'][i], p['y'][i+1]) < y0 and max(p['y'][i], p['y'][i+1]) > y0:
            #print polygon_ordered_list['x'][i+1], x(polygon_ordered_list['y'][i+1])
            if x(y0) > x0:
                point_count_on_the_right += 1
    if point_count_on_the_right % 2 == 0:
        return False
    else:
        return True

def get_density_area_polygon():

    # 2) Show crude plot, so user can click on desired density region
    p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VY].argsort()]
    # Limits with a little margin
    margin = 0.05
    xl = p.max_x - p.min_x; yl = p.max_y - p.min_y
    x0 = p.min_x - margin*xl; xf = p.max_x + margin*xl;
    y0 = p.min_y - margin*yl; yf = p.max_y + margin*yl;
    fig = plt.figure()
    plt.axis([x0, xf, y0, yf])
    plt.axes().set_aspect('equal')
    viewParticles = p.current_matrix[p.current_matrix[:, p.T] == 1]
    sc = plt.scatter(viewParticles[:, p.X], viewParticles[:, p.Y], s=10, linewidth=0)
    viewWall = p.current_matrix[p.current_matrix[:, p.T] == 0]
    sc = plt.scatter(viewWall[:, p.X], viewWall[:, p.Y], s=10, linewidth=0, c = 'red')
    cid = fig.canvas.mpl_connect('button_press_event', plot_onclick)
    cid = fig.canvas.mpl_connect('key_release_event', plot_keyrelease)
    plt.show()

    print "Polygon vertices:"
    print density_polygon_list

def post_analysis_density(import_path):
    # 1) Load state
    saved_states.load_state(import_path)

    # 2) Get the area on where to create the density mesh
    get_density_area_polygon()

    # 3) Create mesh having cells of 0.1 diam x 0.1 diam in the region where we have powder particles (!excluding wall particles!)
    x_min = min(density_polygon_list['x'])
    x_max = max(density_polygon_list['x'])
    y_min = min(density_polygon_list['y'])
    y_max = max(density_polygon_list['y'])

    cell_width = cell_height = 0.1 * p.RADIUS * 2.0
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
            if point_belongs_to_polygon(cell_x, cell_y, density_polygon_list):
                distances = np.square(p.current_matrix[:, p.X:p.Y+1] - np.matrix([cell_x, cell_y]))
                p.current_matrix[:, p.VX] = np.sqrt(distances[:, 0] + distances[:, 1]).transpose()
                # Now we sort by distance (here called VX)
                p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VX].argsort()]
                # The first particle after sorting will have an increment on the number of owned cells
                p.current_matrix[0,p.VY] = p.current_matrix[0,p.VY] + 1

        # Show progress of calculation
        print float(round(float(i)/float(num_x-1)*1000))/10.0, "%"

    # 4) Plot intensity
    p.current_matrix[:,:] = p.current_matrix[p.current_matrix[:,p.VY].argsort()]
    # Limits with a little margin
    margin = 0.05
    xl = x_max - x_min; yl = y_max - y_min
    x0 = x_min - margin*xl; xf = x_max + margin*xl;
    y0 = y_min - margin*yl; yf = y_max + margin*yl;
    plt.axis([x0, xf, y0, yf])
    plt.axes().set_aspect('equal')
    # First filter to eliminate some border effects and division by zero errors
    view = np.logical_and(p.current_matrix[:, p.T] == 1, p.current_matrix[:, p.VY] >= 1)
    dens = (np.pi*(p.RADIUS)**2.0/(p.current_matrix[view, p.VY]*cell_width*cell_height))**0.5
    # Second filter to eliminate the remainder of border effects
    view2 = np.logical_and(dens >= 0.0, dens <= 1.00)
    dens = dens[view2]
    sc = plt.scatter(p.current_matrix[view, p.X][view2], p.current_matrix[view, p.Y][view2], c=dens, s=50, linewidth=0, cmap=plt.cm.rainbow)
    plt.colorbar(sc)
    plt.draw()
    datestr = datetime.datetime.now().isoformat()
    datestr = datestr[0:19].replace(':', '').replace('-', '')
    plt.savefig(os.path.join(import_path, datestr + '_discrete_density_points.pdf'))
    plt.show()

    density = np.zeros((num_x*num_y, 3))
    # 5) Plot Gaussian weighting
    for i in range (0, num_x):
        for j in range(0, num_y):
            cell_x = x_min + i*cell_width
            cell_y = y_min + j*cell_height
            if point_belongs_to_polygon(cell_x, cell_y, density_polygon_list):
                distances = np.square(p.current_matrix[view, p.X:p.Y+1] - np.matrix([cell_x, cell_y]))
                p.current_matrix[view, p.VX] = (distances[:, 0] + distances[:, 1]).transpose()
                # For sigma/d <= 1/4 density variations are resolved at grain level while larger values of sigma/d yield smoother maps(C. Bierwisch et al. / Powder Technology 196 (2009) p. 170)
                sigma = (p.RADIUS*2.0)*1.0/4.0 * 2.0
                dens_gauss = np.sum(dens * np.exp(-p.current_matrix[view, p.VX][view2]/(2.0*sigma**2.0))) / np.sum(np.exp(-p.current_matrix[view, p.VX][view2]/(2.0*sigma**2.0)))
                density[i*num_y + j, :] = np.array((cell_x, cell_y, dens_gauss))
            else:
                density[i*num_y + j, :] = np.array((cell_x, cell_y, None)) # Outside region should not affect the colorbar

        # Show progress of calculation
        print float(round(float(i)/float(num_x-1)*1000))/10.0, "%"

    plt.clf()
    plt.axis([x0, xf, y0, yf])
    plt.axes().set_aspect('equal')
    sc = plt.scatter(density[:,0], density[:,1], c=density[:,2], s=10, linewidth=0, cmap=plt.cm.rainbow)
    plt.colorbar(sc)
    plt.savefig(os.path.join(import_path, datestr + '_gaussian_density_map.png'))
    plt.show()
