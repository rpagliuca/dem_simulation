import numpy as np
import parameters as p

def init():
    # Shoe dimensions
    p.SH = 5.E-2 # shoe height (in m)
    # Multiplicator for random generation
    p.SH_MULTIPLICATOR = 2.
    #L = 4.E-2 # dye length (in m)
    p.SL = p.SH # show lengths (in m)

    # Dye dimensions
    p.DH = 2.E-2
    p.DL = p.DH

    # Number of particles needed to represent the dye (they overlap a little bit ~ 1.6 instead of 2)
    p.NUMBER_PARTICLES_BOTTOM_WALL = np.ceil(p.SL/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_SIDE_WALL = np.ceil(p.SH/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_DYE_BOTTOM_WALL = np.ceil(p.DL/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_DYE_SIDE_WALL = np.ceil(p.DH/(p.RADIUS*1.6))

def calculateN():
    return int(2*p.NUMBER_PARTICLES_BOTTOM_WALL + 2*p.NUMBER_PARTICLES_SIDE_WALL + p.NUMBER_PARTICLES_DYE_BOTTOM_WALL + 2*p.NUMBER_PARTICLES_DYE_SIDE_WALL)

def draw(current_matrix):

    # Table Bottom wall 1
    offset = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.Y] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.X] = np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.WT] = 1 # static wall

    # Table Bottom Wall 2
    offset = offset+p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.Y] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.X] = p.SL + p.DL + np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_BOTTOM_WALL, p.WT] = 1 # static wall

    # Shoe Left wall
    offset = offset+p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.X] = 0
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.WT] = 2 # movable wall

    # Shoe Right wall
    offset = offset+p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.X] = p.SL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_SIDE_WALL, p.WT] = 2 # movable wall

    # Dye Bottom wall
    offset = offset+p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.Y] = -p.DH
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.X] = p.SL + np.linspace(0, p.DL, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL, p.WT] = 1 # static

    # Dye Left Wall
    offset = offset+p.NUMBER_PARTICLES_DYE_BOTTOM_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.X] = p.SL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.WT] = 1 # static

    # Dye Right Wall
    offset = offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.X] = p.SL + p.DL
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.M] = p.WALL_MASS
    current_matrix[offset:offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL, p.WT] = 1 # static

    # Flag particles as type 0 (wall particle)
    offset = offset+p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[0:offset, p.T] = 0

    return current_matrix
