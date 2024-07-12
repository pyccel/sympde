def assemble_matrix_wti60kr7(global_test_basis_v2_1 : "float64[:,:,:,:]", global_test_basis_v2_2 : "float64[:,:,:,:]", global_trial_basis_u2_1 : "float64[:,:,:,:]", global_trial_basis_u2_2 : "float64[:,:,:,:]", global_span_v2_1 : "int64[:]", global_span_v2_2 : "int64[:]", global_x1 : "float64[:,:]", global_x2 : "float64[:,:]", test_v2_p1 : "int64", test_v2_p2 : "int64", trial_u2_p1 : "int64", trial_u2_p2 : "int64", n_element_1 : "int64", n_element_2 : "int64", k1 : "int64", k2 : "int64", pad1 : "int64", pad2 : "int64", g_mat_u2_v2_wti60kr7 : "float64[:,:,:,:]"):

    from numpy import array, zeros, zeros_like, floor
    from math import cos, sqrt, sin, pi
    local_x1 = zeros_like(global_x1[0,:])
    local_x2 = zeros_like(global_x2[0,:])
    
    l_mat_u2_v2_wti60kr7 = zeros((3, 3, 5, 5), dtype='float64')
    for i_element_1 in range(0, n_element_1, 1):
        local_x1[:] = global_x1[i_element_1,:]
        span_v2_1 = global_span_v2_1[i_element_1]
        for i_element_2 in range(0, n_element_2, 1):
            local_x2[:] = global_x2[i_element_2,:]
            span_v2_2 = global_span_v2_2[i_element_2]
            l_mat_u2_v2_wti60kr7[:,:,:,:] = 0.0
            for i_quad_1 in range(0, 4, 1):
                x1 = local_x1[i_quad_1]
                for i_quad_2 in range(0, 4, 1):
                    x2 = local_x2[i_quad_2]
                    for i_basis_1 in range(0, 3, 1):
                        for i_basis_2 in range(0, 3, 1):
                            for j_basis_1 in range(0, 3, 1):
                                v2_1 = global_test_basis_v2_1[i_element_1,i_basis_1,0,i_quad_1]
                                v2_1_x1 = global_test_basis_v2_1[i_element_1,i_basis_1,1,i_quad_1]
                                u2_1 = global_trial_basis_u2_1[i_element_1,j_basis_1,0,i_quad_1]
                                u2_1_x1 = global_trial_basis_u2_1[i_element_1,j_basis_1,1,i_quad_1]
                                for j_basis_2 in range(0, 3, 1):
                                    v2_2 = global_test_basis_v2_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v2_2_x2 = global_test_basis_v2_2[i_element_2,i_basis_2,1,i_quad_2]
                                    u2_2 = global_trial_basis_u2_2[i_element_2,j_basis_2,0,i_quad_2]
                                    u2_2_x2 = global_trial_basis_u2_2[i_element_2,j_basis_2,1,i_quad_2]
                                    v2 = v2_1*v2_2
                                    v2_x2 = v2_1*v2_2_x2
                                    v2_x1 = v2_1_x1*v2_2
                                    u2 = u2_1*u2_2
                                    u2_x2 = u2_1*u2_2_x2
                                    u2_x1 = u2_1_x1*u2_2
                                    temp_v2_u2_0 = 2*pi
                                    temp_v2_u2_1 = temp_v2_u2_0*x1
                                    temp_v2_u2_2 = temp_v2_u2_0*x2
                                    temp_v2_u2_3 = (0.25*sin(temp_v2_u2_1)*cos(temp_v2_u2_2) + 0.25*sin(temp_v2_u2_2)*cos(temp_v2_u2_1) + 1.0)**2
                                    contribution_v2_u2_wti60kr7 = 0.25*u2*v2/sqrt(temp_v2_u2_3)
                                    l_mat_u2_v2_wti60kr7[i_basis_1,i_basis_2,2 - i_basis_1 + j_basis_1,2 - i_basis_2 + j_basis_2] += contribution_v2_u2_wti60kr7
                                
                            
                        
                    
                
            
            g_mat_u2_v2_wti60kr7[pad1 + span_v2_1 - test_v2_p1:1 + pad1 + span_v2_1,pad2 + span_v2_2 - test_v2_p2:1 + pad2 + span_v2_2,:,:] += l_mat_u2_v2_wti60kr7[:,:,:,:]
        
    
    return
