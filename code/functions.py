from numba import jit
from parameters import *

# Source: http://stackoverflow.com/a/29799815/1501575
# Pre-compiled function to find first element of array greater than
@jit(nopython=True)
def find_first_item_greater_than(vec, value):
    for i in xrange(len(vec)):
        if vec[i] > value:
            return i
    return -1

def apply_forces(current_matrix):
    # Now we iterate over every particle, only accounting other particles which y_i - y_j <= 2*RADIUS
    # Last particle shouldn't interact with any other. It has already interacted with the previous ones.
    for i in range (0, N-1):

        # np.argmax returns the minimum index wich satisfies some arbitrary condition, but it is way slower than using numba pre-compiled functions
        #max_neighbour_offset = np.argmax(current_matrix[i+1:, Y] > current_matrix[i, Y] + 2*RADIUS)
        max_neighbour_offset = find_first_item_greater_than(current_matrix[i+1:, Y], current_matrix[i, Y] + 2*RADIUS)

        if max_neighbour_offset == -1:
            max_neighbour_offset = current_matrix[i+1:,Y].size

        max_neighbour_index = i + max_neighbour_offset

        # Particles i through max_neighbour_index are the ones which may interact (based solely on Y distance)
        row_of_interest = current_matrix[i+1:max_neighbour_index+1, :]

        # Now we remove those particles which Y position is greater than radius
        possible_interactions = np.less_equal(row_of_interest[:,X], current_matrix[i, X] + 2*RADIUS)
        cell_of_interest = row_of_interest[possible_interactions, :]

        #contact_forces = 0

        # If there is any possible neighbour
        if cell_of_interest.size > 0:

            # Get distance of possible neighbours
            distances = np.sqrt(np.square(cell_of_interest[:,X] - current_matrix[i,X]) + np.square(cell_of_interest[:,Y] - current_matrix[i, Y]))

            # Generate unitary vector
            radial_unitary_vector = ((cell_of_interest[:, X:Y+1] - current_matrix[i, X:Y+1]).transpose() / (distances + 1.E-20)).transpose()

            # Discard distances greater than 2*RADIUS
            deformations = (2*RADIUS - distances).clip(min=0)

            # Add contact forces

            # Force 1 => Repulsion force
            #contact_forces = (2./3. * E_TILDE * EFFECTIVE_RADIUS**0.5 * deformations**(3./2.) * radial_unitary_vector.transpose()).transpose()
            #contact_forces = (2.*7.48E-7 * 2.0E6 * deformations**(3./2.) * radial_unitary_vector.transpose()).transpose()
            contact_forces = (2.*7.48E-7 * 5.0E-2 / (distances + 1.E-20)**2 * deformations * radial_unitary_vector.transpose()).transpose()

            # Gravity applies force of approx. 9.81*7.63E-8 => 7.48E-7 NEWTON
            # Repulsion should apply more or less the same ammount (maybe double)

            # Force 2 => Attraction (cohesion)
            # Non-existent for single spheres

            # Force 3 => Viscous dissipation
            #contact_forces += - (GAMMA_R * RADIUS * np.sqrt(EFFECTIVE_RADIUS * deformations).transpose() * np.einsum( 'ij, ij->i', (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()
            contact_forces += - 2.0E-5 * (GAMMA_R * RADIUS * np.sqrt(EFFECTIVE_RADIUS * deformations).transpose() * np.einsum( 'ij, ij->i', (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1]) , radial_unitary_vector ) * radial_unitary_vector.transpose()).transpose()

            # Force 4 => Friction force
            # To be modelled
            # I'm not modelling as Cundall, but as Haff and Werner
            relative_velocities = (cell_of_interest[:, VX:VY+1] - current_matrix[i, VX:VY+1])
            tangent_relative_velocities = relative_velocities - (relative_velocities.transpose() - np.einsum('ij, ij->i', relative_velocities, radial_unitary_vector)).transpose() * radial_unitary_vector
            tangent_relative_velocities_module = np.sqrt(tangent_relative_velocities[:, 0]**2 + tangent_relative_velocities[:, 1]**2) + 1.E-20
            tangent_unitary_vector = (tangent_relative_velocities.transpose() / tangent_relative_velocities_module).transpose()
            GBPM_GAMMA = 1.E-6
            #contact_forces += - GBPM_GAMMA * (tangent_relative_velocities.transpose() * deformations).transpose()

            # Write contact_forces to row_of_interest view based on possible_interactions items (indeces)
            row_of_interest[possible_interactions,FX:FY+1] += contact_forces 

            # Add forces and apply its negative sum to current particle (Newton's Second Law of motion)
            current_matrix[i, FX:FY+1] += -np.einsum('ij->j', contact_forces) 
            
    # Apply gravity and Stoke's air drag (except for wall particles by multiplying by column T)
    # Force 5 => Gravity
    # Force 6 => Stoke's air drag
    #current_matrix[:, FX:FY+1] += ((-6 * PI * MU_A * RADIUS * current_matrix[:, VX:VY+1] + (np.outer(current_matrix[:, M]*G, np.array((0, -1))))).transpose() * current_matrix[:, T]).transpose()
    current_matrix[:, FX:FY+1] += (((np.outer(current_matrix[:, M]*G, np.array((0, -1))))).transpose() * current_matrix[:, T]).transpose()

    return current_matrix

