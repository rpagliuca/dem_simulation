# -*- coding: utf-8 -*-

import functions
import numpy as np
import parameters as p # Load all global variables from parameters

# This function checks for particles which are already contacting before the beggining
# of the simulation (due to random position conflict), and removes them
def init_overlap_fix(current_matrix):

    print 'Beggining overlap fix...'

    current_matrix = current_matrix[current_matrix[:,p.Y].argsort()]
    to_be_removed = []

    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= 2*RADIUS
    # Last particle shouldn't interact with any other. It has already interacted with the previous ones.

    for i in range (0, p.N-1):

        max_neighbour_offset = functions.find_first_item_greater_than(current_matrix[i+1:, p.Y], current_matrix[i, p.Y] + 2*p.RADIUS)

        if max_neighbour_offset == -1:
            max_neighbour_offset = current_matrix[i+1:,p.Y].size

        max_neighbour_index = i + max_neighbour_offset

        # Particles i through max_neighbour_index are the ones which may interact (based solely on Y distance)
        row_of_interest = current_matrix[i+1:max_neighbour_index+1, :]

        # Get distance of possible neighbours
        distances = np.sqrt(np.square(row_of_interest[:,p.X] - current_matrix[i,p.X]) + np.square(row_of_interest[:,p.Y] - current_matrix[i, p.Y]))

        # Clip distances to 2*RADIUS, so every touching particle will be set value 0 and stored on clipped_distances
        clipped_distances = (distances - 2*p.RADIUS).clip(min=0)
        touching_particles = row_of_interest[row_of_interest[clipped_distances == 0, p.T] == 1]

        if touching_particles.size > 0:
            # Flag this column as "to be removed" if it is a moving particle (not a wall) and has touching neighbours
            if current_matrix[i, p.T] == 1:
                to_be_removed.append(i)

    current_matrix = np.delete(current_matrix, to_be_removed, 0)
    p.N = len(current_matrix) # update the global variable N, which holds the total of particles

    print ""
    print 'A total of ' + str(len(to_be_removed)) + ' overlapping particles were removed from simulation.'
    print "Number of simulated particles (including wall): " + str(p.N)
    print ""
    print 'Finished overlap fix...'

    return current_matrix
