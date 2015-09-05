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

    # Number of particles needed to represent the dye (they overlap a little bit ~ 1.6 instead of 2)
    p.NUMBER_PARTICLES_BOTTOM_WALL = np.ceil(p.SL/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_SIDE_WALL = np.ceil(p.SH/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_DYE_BOTTOM_WALL = np.ceil(p.DL/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_DYE_SIDE_WALL = np.ceil(p.DH/(p.RADIUS*1.6))
    p.NUMBER_PARTICLES_FLOOR = p.NUMBER_PARTICLES_BOTTOM_WALL*2 + p.NUMBER_PARTICLES_DYE_BOTTOM_WALL

def calculateN():
    return int(2*p.NUMBER_PARTICLES_BOTTOM_WALL + 2*p.NUMBER_PARTICLES_SIDE_WALL + 2*p.NUMBER_PARTICLES_DYE_BOTTOM_WALL + p.NUMBER_PARTICLES_FLOOR + 2*p.NUMBER_PARTICLES_DYE_SIDE_WALL)

def draw(current_matrix):

    # Table Bottom wall 1
    initPos = 0
    endPos = initPos + p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[initPos:endPos, p.Y] = 0
    current_matrix[initPos:endPos, p.X] = np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static wall

    # Table Bottom Wall 2
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_BOTTOM_WALL
    current_matrix[initPos:endPos, p.Y] = 0
    current_matrix[initPos:endPos, p.X] = p.SL + p.DL + np.linspace(0, p.SL, num=p.NUMBER_PARTICLES_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static wall

    # Shoe Left wall
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[initPos:endPos, p.X] = 0
    current_matrix[initPos:endPos, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 2 # movable wall

    # Shoe Right wall
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_SIDE_WALL
    current_matrix[initPos:endPos, p.X] = p.SL
    current_matrix[initPos:endPos, p.Y] = np.linspace(0, p.SH, num=p.NUMBER_PARTICLES_SIDE_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 2 # movable wall

    # Dye Bottom wall (funnel 1)
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_DYE_BOTTOM_WALL
    current_matrix[initPos:endPos, p.Y] = -p.DH - np.linspace(0, p.DL/3, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.X] = p.SL + np.linspace(0, p.DL*0.4, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static

    # Dye Bottom wall (funnel 2)
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_DYE_BOTTOM_WALL
    current_matrix[initPos:endPos, p.Y] = -p.DH - p.DL/3 + np.linspace(0, p.DL/3, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.X] = p.SL + p.DL*0.6 + np.linspace(0, p.DL*0.4, num=p.NUMBER_PARTICLES_DYE_BOTTOM_WALL)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static

    # Floor
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_FLOOR
    current_matrix[initPos:endPos, p.Y] = -p.DH * 3.0
    current_matrix[initPos:endPos, p.X] = np.linspace(0, p.DL + 2*p.SL, num=p.NUMBER_PARTICLES_FLOOR)
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static

    # Dye Left Wall
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[initPos:endPos, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[initPos:endPos, p.X] = p.SL
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static

    # Dye Right Wall
    initPos = endPos
    endPos = initPos + p.NUMBER_PARTICLES_DYE_SIDE_WALL
    current_matrix[initPos:endPos, p.Y] = 0 - np.linspace(0, p.DH, num=p.NUMBER_PARTICLES_DYE_SIDE_WALL)
    current_matrix[initPos:endPos, p.X] = p.SL + p.DL
    current_matrix[initPos:endPos, p.M] = p.WALL_MASS
    current_matrix[initPos:endPos, p.WT] = 1 # static

    # Flag particles as type 0 (wall particle)
    current_matrix[0:endPos, p.T] = 0

    return current_matrix
