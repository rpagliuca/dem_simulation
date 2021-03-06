import numpy as np
import parameters as p

def init():
    # Shoe dimensions
    p.SH = 5.E-2 # shoe height (in m)
    # Multiplicator for random generation
    p.SH_MULTIPLICATOR = 2.
    p.SL = p.SH # shoe lengths (in m)

    # Dye dimensions
    p.DH = 2.E-2
    p.DL = p.DH

def calculateN():
    global endPos
    draw([], True) # Call draw() with calculateN = True as parameter
    return int(endPos)

def draw(current_matrix, calculateN = False):

    global endPos
    endPos = 0

    # Table Bottom wall 1
    draw_line(0, 0, p.SL, 0, current_matrix,1, calculateN)
    # Table Bottom Wall 2
    draw_line(p.SL + p.DL, 0, p.SL + p.DL + p.SL, 0, current_matrix, 1, calculateN)
    # Shoe Left wall (movable wall)
    draw_line(0, 0, 0, p.SH, current_matrix, 2, calculateN) # wall_type -> movable_wall
    # Shoe Right wall (movable wall)
    draw_line(p.SL, 0, p.SL, p.SH, current_matrix, 2, calculateN) # wall_type -> movable_wall
    # Dye Left Wall
    draw_line(p.SL, 0, p.SL, -p.DH, current_matrix, 1, calculateN)
    # Dye Right Wall
    draw_line(p.SL + p.DL, 0, p.SL + p.DL, -p.DH, current_matrix, 1, calculateN)
    # Funnel left slope
    draw_line(p.SL, -p.DH, p.SL + p.DL*0.4, -p.DH - p.DL/3.0, current_matrix, 1, calculateN)
    # Funnel left straight tip
    draw_line(p.SL + p.DL*0.4, -p.DH - p.DL/3.0, p.SL + p.DL*0.4, -p.DH - p.DL/3.0 - p.DL*0.2, current_matrix, 1, calculateN)
    # Funnel right slope
    draw_line(p.SL + p.DL*0.6, -p.DH - p.DL/3.0, p.SL + p.DL*0.6 + p.DL*0.4, -p.DH - p.DL/3.0 + p.DL/3.0, current_matrix, 1, calculateN) 
    # Funnel right straight tip
    draw_line(p.SL + p.DL*0.6, -p.DH - p.DL/3.0, p.SL + p.DL*0.6, -p.DH - p.DL/3.0 - p.DL*0.2, current_matrix, 1, calculateN)
    # Floor Bottom
    x_center = (p.DL + 2*p.SL)/2.0
    draw_line(x_center - p.DL, -p.DH * 2.5, x_center + p.DL, -p.DH * 2.5, current_matrix, 1, calculateN)
    # Floor Left Wall
    draw_line(x_center - p.DL, -p.DH * 2.5, x_center - p.DL, -p.DH * 2.0, current_matrix, 1, calculateN)
    # Floor Right Wall
    draw_line(x_center + p.DL, -p.DH * 2.5, x_center + p.DL, -p.DH * 2.0, current_matrix, 1, calculateN)

    if not calculateN:
        # Common properties for wall particles
        current_matrix[0:endPos, p.T] = 0 # Type 0 -> Wall particle
        current_matrix[0:endPos, p.M] = p.WALL_MASS

    return current_matrix

def draw_line(x0, y0, xf, yf, current_matrix, wall_type, calculateN):

    global endPos

    # Calculate length of line and define number of particles needed
    number_particles = np.ceil((((x0-xf)**2 + (y0-yf)**2)**0.5)/(p.RADIUS*1.6))
    initPos = endPos
    endPos = initPos + number_particles

    if not calculateN:

        # Draw particles as wall particles
        current_matrix[initPos:endPos, p.Y] = y0 + np.linspace(0, yf-y0, num=number_particles)
        current_matrix[initPos:endPos, p.X] = x0 + np.linspace(0, xf-x0, num=number_particles)
        current_matrix[initPos:endPos, p.WT] = wall_type

        # Update plot limits
        if x0 < p.min_x:
            p.min_x = x0
        if xf < p.min_x:
            p.min_x = xf
        if x0 > p.max_x:
            p.max_x = x0
        if xf > p.max_x:
            p.max_x = xf
        if y0 < p.min_y:
            p.min_y = y0
        if yf < p.min_y:
            p.min_y = yf
        if y0 > p.max_y:
            p.max_y = y0
        if yf > p.max_y:
            p.max_y = yf
