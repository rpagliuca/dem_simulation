#radial_unitary_vector = np.array([[np.sqrt(0.5), np.sqrt(0.5)], [np.sqrt(0.5), np.sqrt(0.5)], [np.sqrt(0.5), np.sqrt(0.5)]])
radial_unitary_vector = np.array([[0, 1], [0, 1], [0, 1]])
current_matrix[i, VX:VY+1] = np.array((2., 4.))
cell_of_interest = np.array([[4., 6.], [1., 1.], [2., 2.]])
#print current_matrix[i, VX:VY+1]
#print
#print cell_of_interest
#print
vel_rel = cell_of_interest[:, 0:2] - current_matrix[i, VX:VY+1]
print vel_rel  
print
print radial_unitary_vector
print
projection = np.einsum( 'ij, ij->i', vel_rel, radial_unitary_vector)
print projection
print
print (projection * radial_unitary_vector.transpose()).transpose()

exit()

