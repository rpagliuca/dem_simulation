# -*- coding: utf-8 -*-

import functions
from parameters import * # Load all global variables from parameters

def init_overlap_fix(current_matrix):

    print 'Beggining overlap fix...'

    # First, we sort by y position to optimize contact forces
    current_matrix = current_matrix[current_matrix[:,Y].argsort()]

    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= 2*RADIUS
    # Last particle shouldn't interact with any other. It has already interacted with the previous ones.

    for i in range (0, N-1):

        max_neighbour_offset = functions.find_first_item_greater_than(current_matrix[i+1:, Y], current_matrix[i, Y] + 2*RADIUS)

        if max_neighbour_offset == -1:
            max_neighbour_offset = current_matrix[i+1:,Y].size

        max_neighbour_index = i + max_neighbour_offset

        # Particles i through max_neighbour_index are the ones which may interact (based solely on Y distance)
        row_of_interest = current_matrix[i+1:max_neighbour_index+1, :]

        # Get distance of possible neighbours
        distances = np.sqrt(np.square(row_of_interest[:,X] - current_matrix[i,X]) + np.square(row_of_interest[:,Y] - current_matrix[i, Y]))

        # Clip distances to 2*RADIUS, so every touching particle will be set value 0 and stored on clipped_distances
        clipped_distances = (distances - 2*RADIUS).clip(min=0)
        touching_particles = row_of_interest[row_of_interest[clipped_distances == 0, T] == 1]

        if touching_particles.size > 0:
            # Zero everything from this column if it is a particle and has touching neighbours
            if current_matrix[i, T] == 1:
                current_matrix[i, :] = 0

    print 'Finished overlap fix...'

    return current_matrix
