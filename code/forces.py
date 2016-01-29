# -*- coding: utf-8 -*-

import functions
import parameters as p # Load all global variables from parameters
import numpy as np
import joblib
import os

def apply_forces(current_matrix, parallel):

    # Using multicores, we apply forces that depend on the neighbours of the particle (we must loop through them all)
    n_split = int(np.floor(p.N/p.number_of_cores))
    parallel( # Multi core
        joblib.delayed(loop_forces)(current_matrix, i*n_split, (i+1)*n_split-1) for i in range(0, p.number_of_cores))

    # Apply forces that do not depend on neighbours
    # Gravity and Stoke's air drag (except for wall particles -- by multiplying by column T we have null force for wall particles)
    # Force 5 => Gravity
    # and
    # Force 6 => Stoke's air drag
    current_matrix[:, p.FX:p.FY+1] += ((-6 * p.PI * p.MU_A * p.RADIUS * current_matrix[:, p.VX:p.VY+1] + (np.outer(current_matrix[:, p.M]*p.G, np.array((0, -1))))).transpose() * current_matrix[:, p.T]).transpose()
    current_matrix[:, p.FX:p.FY+1] += (((np.outer(current_matrix[:, p.M]*p.G, np.array((0, -1))))).transpose() * current_matrix[:, p.T]).transpose()
    return current_matrix

def loop_forces(current_matrix, start_index = False, end_index = False):

    # This function MUST receive a previously sorted matrix based on Y position; failing to do so, it will return the wrong results

    # Default indices
    if not start_index:
        start_index = 0
    if not end_index:
        end_index = N-1

    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= 2*RADIUS
    # Last particle shouldn't interact with any other. It has already interacted with the previous ones.
    for i in range (start_index, end_index):

        # np.argmax returns the minimum index wich satisfies some arbitrary condition, but it is way slower than using numba pre-compiled functions
        #max_neighbour_offset = np.argmax(current_matrix[i+1:, Y] > current_matrix[i, Y] + 2*RADIUS)
        max_neighbour_offset = functions.find_first_item_greater_than(current_matrix[i+1:, p.Y], current_matrix[i, p.Y] + 2*p.RADIUS)

        if max_neighbour_offset == -1:
            max_neighbour_offset = current_matrix[i+1:,p.Y].size

        max_neighbour_index = i + max_neighbour_offset

        # Particles i through max_neighbour_index are the ones which may interact (based solely on Y distance)
        row_of_interest = current_matrix[i+1:max_neighbour_index+1, :]

        # Now we remove those particles which X position is greater than radius
        possible_interactions = np.less_equal(row_of_interest[:,p.X], current_matrix[i, p.X] + 2*p.RADIUS)
        cell_of_interest = row_of_interest[possible_interactions, :]

        # If there is any possible neighbour
        if cell_of_interest.size > 0:

            # Get distances of possible neighbours
            distances = np.sqrt(np.square(cell_of_interest[:,p.X] - current_matrix[i,p.X]) + np.square(cell_of_interest[:,p.Y] - current_matrix[i, p.Y]))

            # Define the unitary vector
            radial_unitary_vector = ((cell_of_interest[:, p.X:p.Y+1] - current_matrix[i, p.X:p.Y+1]).transpose() / (distances + 1.E-20)).transpose()

            # Discard distances greater than 2*RADIUS
            deformations = (2*p.RADIUS - distances).clip(min=0)

            # Add contact forces

            # Force 1 => Repulsion force
            contact_forces = (2./3. * p.E_TILDE * p.EFFECTIVE_RADIUS**0.5 * deformations**(3./2.) * radial_unitary_vector.transpose()).transpose()

            # Force 2 => Attraction (cohesion)
            # Non-existent for single spheres

            # Force 3 => Viscous dissipation
            contact_forces += - (p.GAMMA_R * p.RADIUS * np.sqrt(p.EFFECTIVE_RADIUS * deformations).transpose() * np.einsum( 'ij, ij->i', (cell_of_interest[:, p.VX:p.VY+1] - current_matrix[i, p.VX:p.VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()

            # Force 4 => Friction force
            # I'm not modelling this force as Cundall, but as Haff and Werner
            relative_velocities = (cell_of_interest[:, p.VX:p.VY+1] - current_matrix[i, p.VX:p.VY+1])
            tangent_relative_velocities = relative_velocities - (relative_velocities.transpose() - np.einsum('ij, ij->i', relative_velocities, radial_unitary_vector)).transpose() * radial_unitary_vector
            tangent_relative_velocities_module = np.sqrt(tangent_relative_velocities[:, 0]**2 + tangent_relative_velocities[:, 1]**2) + 1.E-20
            tangent_unitary_vector = (tangent_relative_velocities.transpose() / tangent_relative_velocities_module).transpose()
            contact_forces += - p.GBPM_GAMMA * (tangent_relative_velocities.transpose() * deformations).transpose()
            
            # Write contact_forces to row_of_interest view based on possible_interactions items (indices)
            row_of_interest[possible_interactions,p.FX:p.FY+1] += contact_forces 

            # Add forces and apply its negative sum to current particle (Newton's Third Law of motion)
            current_matrix[i, p.FX:p.FY+1] += -np.einsum('ij->j', contact_forces) 
