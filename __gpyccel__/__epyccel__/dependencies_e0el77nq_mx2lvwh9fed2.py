def assemble_matrix_e0el77nq(global_test_basis_v1_0_1 : "float64[:,:,:,:]", global_test_basis_v1_0_2 : "float64[:,:,:,:]", global_test_basis_v1_1_1 : "float64[:,:,:,:]", global_test_basis_v1_1_2 : "float64[:,:,:,:]", global_trial_basis_u1_0_1 : "float64[:,:,:,:]", global_trial_basis_u1_0_2 : "float64[:,:,:,:]", global_trial_basis_u1_1_1 : "float64[:,:,:,:]", global_trial_basis_u1_1_2 : "float64[:,:,:,:]", global_span_v1_0_1 : "int64[:]", global_span_v1_0_2 : "int64[:]", global_span_v1_1_1 : "int64[:]", global_span_v1_1_2 : "int64[:]", global_x1 : "float64[:,:]", global_x2 : "float64[:,:]", test_v1_0_p1 : "int64", test_v1_0_p2 : "int64", test_v1_1_p1 : "int64", test_v1_1_p2 : "int64", trial_u1_0_p1 : "int64", trial_u1_0_p2 : "int64", trial_u1_1_p1 : "int64", trial_u1_1_p2 : "int64", n_element_1 : "int64", n_element_2 : "int64", k1 : "int64", k2 : "int64", pad1 : "int64", pad2 : "int64", g_mat_u1_0_v1_0_e0el77nq : "float64[:,:,:,:]", g_mat_u1_1_v1_0_e0el77nq : "float64[:,:,:,:]", g_mat_u1_0_v1_1_e0el77nq : "float64[:,:,:,:]", g_mat_u1_1_v1_1_e0el77nq : "float64[:,:,:,:]"):

    from numpy import array, zeros, zeros_like, floor
    from math import cos, sqrt, sin, pi
    local_x1 = zeros_like(global_x1[0,:])
    local_x2 = zeros_like(global_x2[0,:])
    
    l_mat_u1_0_v1_0_e0el77nq = zeros((3, 4, 5, 7), dtype='float64')
    l_mat_u1_0_v1_1_e0el77nq = zeros((4, 3, 7, 7), dtype='float64')
    l_mat_u1_1_v1_0_e0el77nq = zeros((3, 4, 7, 7), dtype='float64')
    l_mat_u1_1_v1_1_e0el77nq = zeros((4, 3, 7, 5), dtype='float64')
    for i_element_1 in range(0, n_element_1, 1):
        local_x1[:] = global_x1[i_element_1,:]
        span_v1_0_1 = global_span_v1_0_1[i_element_1]
        span_v1_1_1 = global_span_v1_1_1[i_element_1]
        for i_element_2 in range(0, n_element_2, 1):
            local_x2[:] = global_x2[i_element_2,:]
            span_v1_0_2 = global_span_v1_0_2[i_element_2]
            span_v1_1_2 = global_span_v1_1_2[i_element_2]
            l_mat_u1_0_v1_0_e0el77nq[:,:,:,:] = 0.0
            l_mat_u1_1_v1_0_e0el77nq[:,:,:,:] = 0.0
            l_mat_u1_0_v1_1_e0el77nq[:,:,:,:] = 0.0
            l_mat_u1_1_v1_1_e0el77nq[:,:,:,:] = 0.0
            for i_quad_1 in range(0, 4, 1):
                x1 = local_x1[i_quad_1]
                for i_quad_2 in range(0, 4, 1):
                    x2 = local_x2[i_quad_2]
                    for i_basis_1 in range(0, 3, 1):
                        for i_basis_2 in range(0, 4, 1):
                            for j_basis_1 in range(0, 3, 1):
                                v1_0_1 = global_test_basis_v1_0_1[i_element_1,i_basis_1,0,i_quad_1]
                                v1_0_1_x1 = global_test_basis_v1_0_1[i_element_1,i_basis_1,1,i_quad_1]
                                u1_0_1 = global_trial_basis_u1_0_1[i_element_1,j_basis_1,0,i_quad_1]
                                u1_0_1_x1 = global_trial_basis_u1_0_1[i_element_1,j_basis_1,1,i_quad_1]
                                for j_basis_2 in range(0, 4, 1):
                                    v1_0_2 = global_test_basis_v1_0_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v1_0_2_x2 = global_test_basis_v1_0_2[i_element_2,i_basis_2,1,i_quad_2]
                                    u1_0_2 = global_trial_basis_u1_0_2[i_element_2,j_basis_2,0,i_quad_2]
                                    u1_0_2_x2 = global_trial_basis_u1_0_2[i_element_2,j_basis_2,1,i_quad_2]
                                    v1_0 = v1_0_1*v1_0_2
                                    v1_0_x2 = v1_0_1*v1_0_2_x2
                                    v1_0_x1 = v1_0_1_x1*v1_0_2
                                    u1_0 = u1_0_1*u1_0_2
                                    u1_0_x2 = u1_0_1*u1_0_2_x2
                                    u1_0_x1 = u1_0_1_x1*u1_0_2
                                    temp_v1_0_u1_0_0 = 2*pi
                                    temp_v1_0_u1_0_1 = temp_v1_0_u1_0_0*x2
                                    temp_v1_0_u1_0_2 = temp_v1_0_u1_0_0*x1
                                    temp_v1_0_u1_0_3 = 0.25*sin(temp_v1_0_u1_0_1)*cos(temp_v1_0_u1_0_2)
                                    temp_v1_0_u1_0_4 = sin(temp_v1_0_u1_0_2)
                                    temp_v1_0_u1_0_5 = cos(temp_v1_0_u1_0_1)
                                    temp_v1_0_u1_0_6 = 0.25*temp_v1_0_u1_0_4*temp_v1_0_u1_0_5
                                    temp_v1_0_u1_0_7 = temp_v1_0_u1_0_6 + 1.0
                                    temp_v1_0_u1_0_8 = u1_0*v1_0/(temp_v1_0_u1_0_3 + temp_v1_0_u1_0_6 + 1)**2
                                    contribution_v1_0_u1_0_e0el77nq = 4.0*(0.015625*temp_v1_0_u1_0_4**2*temp_v1_0_u1_0_5**2*temp_v1_0_u1_0_8 + 0.25*temp_v1_0_u1_0_7**2*temp_v1_0_u1_0_8)*sqrt((temp_v1_0_u1_0_3 + temp_v1_0_u1_0_7)**2)
                                    l_mat_u1_0_v1_0_e0el77nq[i_basis_1,i_basis_2,2 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2] += contribution_v1_0_u1_0_e0el77nq
                                
                            
                        
                    
                    for i_basis_1 in range(0, 3, 1):
                        for i_basis_2 in range(0, 4, 1):
                            for j_basis_1 in range(0, 4, 1):
                                v1_0_1 = global_test_basis_v1_0_1[i_element_1,i_basis_1,0,i_quad_1]
                                v1_0_1_x1 = global_test_basis_v1_0_1[i_element_1,i_basis_1,1,i_quad_1]
                                u1_1_1 = global_trial_basis_u1_1_1[i_element_1,j_basis_1,0,i_quad_1]
                                u1_1_1_x1 = global_trial_basis_u1_1_1[i_element_1,j_basis_1,1,i_quad_1]
                                for j_basis_2 in range(0, 3, 1):
                                    v1_0_2 = global_test_basis_v1_0_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v1_0_2_x2 = global_test_basis_v1_0_2[i_element_2,i_basis_2,1,i_quad_2]
                                    u1_1_2 = global_trial_basis_u1_1_2[i_element_2,j_basis_2,0,i_quad_2]
                                    u1_1_2_x2 = global_trial_basis_u1_1_2[i_element_2,j_basis_2,1,i_quad_2]
                                    v1_0 = v1_0_1*v1_0_2
                                    v1_0_x2 = v1_0_1*v1_0_2_x2
                                    v1_0_x1 = v1_0_1_x1*v1_0_2
                                    u1_1 = u1_1_1*u1_1_2
                                    u1_1_x2 = u1_1_1*u1_1_2_x2
                                    u1_1_x1 = u1_1_1_x1*u1_1_2
                                    temp_v1_0_u1_1_0 = 2*pi
                                    temp_v1_0_u1_1_1 = temp_v1_0_u1_1_0*x1
                                    temp_v1_0_u1_1_2 = temp_v1_0_u1_1_0*x2
                                    temp_v1_0_u1_1_3 = sin(temp_v1_0_u1_1_2)*cos(temp_v1_0_u1_1_1)
                                    temp_v1_0_u1_1_4 = sin(temp_v1_0_u1_1_1)*cos(temp_v1_0_u1_1_2)
                                    temp_v1_0_u1_1_5 = 0.25*temp_v1_0_u1_1_4
                                    temp_v1_0_u1_1_6 = temp_v1_0_u1_1_5 + 1.0
                                    temp_v1_0_u1_1_7 = 0.5*temp_v1_0_u1_1_3
                                    temp_v1_0_u1_1_8 = temp_v1_0_u1_1_7 + 2.0
                                    temp_v1_0_u1_1_9 = u1_1*v1_0/((0.5*temp_v1_0_u1_1_4 + temp_v1_0_u1_1_8)*(1.0*temp_v1_0_u1_1_3 + 1.0*temp_v1_0_u1_1_4 + 4.0))
                                    contribution_v1_0_u1_1_e0el77nq = 4.0*(-temp_v1_0_u1_1_5*temp_v1_0_u1_1_8*temp_v1_0_u1_1_9 - temp_v1_0_u1_1_6*temp_v1_0_u1_1_7*temp_v1_0_u1_1_9)*sqrt((0.25*temp_v1_0_u1_1_3 + temp_v1_0_u1_1_6)**2)
                                    l_mat_u1_1_v1_0_e0el77nq[i_basis_1,i_basis_2,3 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2] += contribution_v1_0_u1_1_e0el77nq
                                
                            
                        
                    
                    for i_basis_1 in range(0, 4, 1):
                        for i_basis_2 in range(0, 3, 1):
                            for j_basis_1 in range(0, 3, 1):
                                v1_1_1 = global_test_basis_v1_1_1[i_element_1,i_basis_1,0,i_quad_1]
                                v1_1_1_x1 = global_test_basis_v1_1_1[i_element_1,i_basis_1,1,i_quad_1]
                                u1_0_1 = global_trial_basis_u1_0_1[i_element_1,j_basis_1,0,i_quad_1]
                                u1_0_1_x1 = global_trial_basis_u1_0_1[i_element_1,j_basis_1,1,i_quad_1]
                                for j_basis_2 in range(0, 4, 1):
                                    v1_1_2 = global_test_basis_v1_1_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v1_1_2_x2 = global_test_basis_v1_1_2[i_element_2,i_basis_2,1,i_quad_2]
                                    u1_0_2 = global_trial_basis_u1_0_2[i_element_2,j_basis_2,0,i_quad_2]
                                    u1_0_2_x2 = global_trial_basis_u1_0_2[i_element_2,j_basis_2,1,i_quad_2]
                                    v1_1 = v1_1_1*v1_1_2
                                    v1_1_x2 = v1_1_1*v1_1_2_x2
                                    v1_1_x1 = v1_1_1_x1*v1_1_2
                                    u1_0 = u1_0_1*u1_0_2
                                    u1_0_x2 = u1_0_1*u1_0_2_x2
                                    u1_0_x1 = u1_0_1_x1*u1_0_2
                                    temp_v1_1_u1_0_0 = 2*pi
                                    temp_v1_1_u1_0_1 = temp_v1_1_u1_0_0*x1
                                    temp_v1_1_u1_0_2 = temp_v1_1_u1_0_0*x2
                                    temp_v1_1_u1_0_3 = sin(temp_v1_1_u1_0_2)*cos(temp_v1_1_u1_0_1)
                                    temp_v1_1_u1_0_4 = sin(temp_v1_1_u1_0_1)*cos(temp_v1_1_u1_0_2)
                                    temp_v1_1_u1_0_5 = 0.25*temp_v1_1_u1_0_4
                                    temp_v1_1_u1_0_6 = temp_v1_1_u1_0_5 + 1.0
                                    temp_v1_1_u1_0_7 = 0.5*temp_v1_1_u1_0_3
                                    temp_v1_1_u1_0_8 = temp_v1_1_u1_0_7 + 2.0
                                    temp_v1_1_u1_0_9 = u1_0*v1_1/((0.5*temp_v1_1_u1_0_4 + temp_v1_1_u1_0_8)*(1.0*temp_v1_1_u1_0_3 + 1.0*temp_v1_1_u1_0_4 + 4.0))
                                    contribution_v1_1_u1_0_e0el77nq = 4.0*(-temp_v1_1_u1_0_5*temp_v1_1_u1_0_8*temp_v1_1_u1_0_9 - temp_v1_1_u1_0_6*temp_v1_1_u1_0_7*temp_v1_1_u1_0_9)*sqrt((0.25*temp_v1_1_u1_0_3 + temp_v1_1_u1_0_6)**2)
                                    l_mat_u1_0_v1_1_e0el77nq[i_basis_1,i_basis_2,3 - i_basis_1 + j_basis_1,3 - i_basis_2 + j_basis_2] += contribution_v1_1_u1_0_e0el77nq
                                
                            
                        
                    
                    for i_basis_1 in range(0, 4, 1):
                        for i_basis_2 in range(0, 3, 1):
                            for j_basis_1 in range(0, 4, 1):
                                v1_1_1 = global_test_basis_v1_1_1[i_element_1,i_basis_1,0,i_quad_1]
                                v1_1_1_x1 = global_test_basis_v1_1_1[i_element_1,i_basis_1,1,i_quad_1]
                                u1_1_1 = global_trial_basis_u1_1_1[i_element_1,j_basis_1,0,i_quad_1]
                                u1_1_1_x1 = global_trial_basis_u1_1_1[i_element_1,j_basis_1,1,i_quad_1]
                                for j_basis_2 in range(0, 3, 1):
                                    v1_1_2 = global_test_basis_v1_1_2[i_element_2,i_basis_2,0,i_quad_2]
                                    v1_1_2_x2 = global_test_basis_v1_1_2[i_element_2,i_basis_2,1,i_quad_2]
                                    u1_1_2 = global_trial_basis_u1_1_2[i_element_2,j_basis_2,0,i_quad_2]
                                    u1_1_2_x2 = global_trial_basis_u1_1_2[i_element_2,j_basis_2,1,i_quad_2]
                                    v1_1 = v1_1_1*v1_1_2
                                    v1_1_x2 = v1_1_1*v1_1_2_x2
                                    v1_1_x1 = v1_1_1_x1*v1_1_2
                                    u1_1 = u1_1_1*u1_1_2
                                    u1_1_x2 = u1_1_1*u1_1_2_x2
                                    u1_1_x1 = u1_1_1_x1*u1_1_2
                                    temp_v1_1_u1_1_0 = 2*pi
                                    temp_v1_1_u1_1_1 = temp_v1_1_u1_1_0*x1
                                    temp_v1_1_u1_1_2 = temp_v1_1_u1_1_0*x2
                                    temp_v1_1_u1_1_3 = 0.25*sin(temp_v1_1_u1_1_1)*cos(temp_v1_1_u1_1_2)
                                    temp_v1_1_u1_1_4 = sin(temp_v1_1_u1_1_2)
                                    temp_v1_1_u1_1_5 = cos(temp_v1_1_u1_1_1)
                                    temp_v1_1_u1_1_6 = 0.25*temp_v1_1_u1_1_4*temp_v1_1_u1_1_5
                                    temp_v1_1_u1_1_7 = temp_v1_1_u1_1_6 + 1
                                    temp_v1_1_u1_1_8 = u1_1*v1_1/(temp_v1_1_u1_1_3 + temp_v1_1_u1_1_7)**2
                                    contribution_v1_1_u1_1_e0el77nq = 4.0*(0.015625*temp_v1_1_u1_1_4**2*temp_v1_1_u1_1_5**2*temp_v1_1_u1_1_8 + 0.25*temp_v1_1_u1_1_7**2*temp_v1_1_u1_1_8)*sqrt((temp_v1_1_u1_1_3 + temp_v1_1_u1_1_6 + 1.0)**2)
                                    l_mat_u1_1_v1_1_e0el77nq[i_basis_1,i_basis_2,3 - i_basis_1 + j_basis_1,2 - i_basis_2 + j_basis_2] += contribution_v1_1_u1_1_e0el77nq
                                
                            
                        
                    
                
            
            g_mat_u1_0_v1_0_e0el77nq[pad1 + span_v1_0_1 - test_v1_0_p1:1 + pad1 + span_v1_0_1,pad2 + span_v1_0_2 - test_v1_0_p2:1 + pad2 + span_v1_0_2,:,:] += l_mat_u1_0_v1_0_e0el77nq[:,:,:,:]
            g_mat_u1_1_v1_0_e0el77nq[pad1 + span_v1_0_1 - test_v1_0_p1:1 + pad1 + span_v1_0_1,pad2 + span_v1_0_2 - test_v1_0_p2:1 + pad2 + span_v1_0_2,:,:] += l_mat_u1_1_v1_0_e0el77nq[:,:,:,:]
            g_mat_u1_0_v1_1_e0el77nq[pad1 + span_v1_1_1 - test_v1_1_p1:1 + pad1 + span_v1_1_1,pad2 + span_v1_1_2 - test_v1_1_p2:1 + pad2 + span_v1_1_2,:,:] += l_mat_u1_0_v1_1_e0el77nq[:,:,:,:]
            g_mat_u1_1_v1_1_e0el77nq[pad1 + span_v1_1_1 - test_v1_1_p1:1 + pad1 + span_v1_1_1,pad2 + span_v1_1_2 - test_v1_1_p2:1 + pad2 + span_v1_1_2,:,:] += l_mat_u1_1_v1_1_e0el77nq[:,:,:,:]
        
    
    return
