def lo_dot_vga8r193(mat00 : "float64[:,:,:,:]", mat01 : "float64[:,:,:,:]", mat10 : "float64[:,:,:,:]", mat11 : "float64[:,:,:,:]", x0 : "float64[:,:]", x1 : "float64[:,:]", out0 : "float64[:,:]", out1 : "float64[:,:]", s00_1 : "int64", s00_2 : "int64", s01_1 : "int64", s01_2 : "int64", s10_1 : "int64", s10_2 : "int64", s11_1 : "int64", s11_2 : "int64", n00_1 : "int64", n00_2 : "int64", n01_1 : "int64", n01_2 : "int64", n10_1 : "int64", n10_2 : "int64", n11_1 : "int64", n11_2 : "int64", ne00_1 : "int64", ne00_2 : "int64", ne01_1 : "int64", ne01_2 : "int64", ne10_1 : "int64", ne10_2 : "int64", ne11_1 : "int64", ne11_2 : "int64"):

    
    for i1 in range(0, n00_1, 1):
        for i2 in range(0, n00_2, 1):
            v00 = 0.0
            for k1 in range(0, 5, 1):
                for k2 in range(0, 7, 1):
                    v00 += mat00[3 + i1,3 + i2,k1,k2]*x0[1 + i1 + k1,i2 + k2]
                
            
            out0[3 + i1,3 + i2] = v00
        
    
    for i1 in range(0, ne00_1, 1):
        for i2 in range(0, n00_2, 1):
            v00 = 0.0
            for k1 in range(0, 4 - i1, 1):
                for k2 in range(0, 7, 1):
                    v00 += x0[1 + i1 + k1 + n00_1,i2 + k2]*mat00[3 + i1 + n00_1,3 + i2,k1,k2]
                
            
            out0[3 + i1 + n00_1,3 + i2] = v00
        
    
    for i1 in range(0, n00_1 + ne00_1, 1):
        for i2 in range(0, ne00_2, 1):
            v00 = 0.0
            for k1 in range(0, 5 - max(0, i1 - n00_1 + 1), 1):
                for k2 in range(0, 6 - i2, 1):
                    v00 += x0[1 + i1 + k1,i2 + k2 + n00_2]*mat00[3 + i1,3 + i2 + n00_2,k1,k2]
                
            
            out0[3 + i1,3 + i2 + n00_2] = v00
        
    
    for i1 in range(0, n11_1, 1):
        for i2 in range(0, n11_2, 1):
            v11 = 0.0
            for k1 in range(0, 7, 1):
                for k2 in range(0, 5, 1):
                    v11 += mat11[3 + i1,3 + i2,k1,k2]*x1[i1 + k1,1 + i2 + k2]
                
            
            out1[3 + i1,3 + i2] = v11
        
    
    for i1 in range(0, ne11_1, 1):
        for i2 in range(0, n11_2, 1):
            v11 = 0.0
            for k1 in range(0, 6 - i1, 1):
                for k2 in range(0, 5, 1):
                    v11 += x1[i1 + k1 + n11_1,1 + i2 + k2]*mat11[3 + i1 + n11_1,3 + i2,k1,k2]
                
            
            out1[3 + i1 + n11_1,3 + i2] = v11
        
    
    for i1 in range(0, n11_1 + ne11_1, 1):
        for i2 in range(0, ne11_2, 1):
            v11 = 0.0
            for k1 in range(0, 7 - max(0, i1 - n11_1 + 1), 1):
                for k2 in range(0, 4 - i2, 1):
                    v11 += x1[i1 + k1,1 + i2 + k2 + n11_2]*mat11[3 + i1,3 + i2 + n11_2,k1,k2]
                
            
            out1[3 + i1,3 + i2 + n11_2] = v11
        
    
    for i1 in range(0, n01_1, 1):
        for i2 in range(0, n01_2, 1):
            v01 = 0.0
            for k1 in range(0, 7, 1):
                for k2 in range(0, 7, 1):
                    v01 += mat01[3 + i1,3 + i2,k1,k2]*x1[i1 + k1,i2 + k2]
                
            
            out0[3 + i1,3 + i2] += v01
        
    
    for i1 in range(0, ne01_1, 1):
        for i2 in range(0, n01_2, 1):
            v01 = 0.0
            for k1 in range(0, 6 - i1, 1):
                for k2 in range(0, 7, 1):
                    v01 += x1[i1 + k1 + n01_1,i2 + k2]*mat01[3 + i1 + n01_1,3 + i2,k1,k2]
                
            
            out0[3 + i1 + n01_1,3 + i2] += v01
        
    
    for i1 in range(0, n01_1 + ne01_1, 1):
        for i2 in range(0, ne01_2, 1):
            v01 = 0.0
            for k1 in range(0, 7 - max(0, i1 - n01_1 + 1), 1):
                for k2 in range(0, 6 - i2, 1):
                    v01 += x1[i1 + k1,i2 + k2 + n01_2]*mat01[3 + i1,3 + i2 + n01_2,k1,k2]
                
            
            out0[3 + i1,3 + i2 + n01_2] += v01
        
    
    for i1 in range(0, n10_1, 1):
        for i2 in range(0, n10_2, 1):
            v10 = 0.0
            for k1 in range(0, 7, 1):
                for k2 in range(0, 7, 1):
                    v10 += mat10[3 + i1,3 + i2,k1,k2]*x0[i1 + k1,i2 + k2]
                
            
            out1[3 + i1,3 + i2] += v10
        
    
    for i1 in range(0, ne10_1, 1):
        for i2 in range(0, n10_2, 1):
            v10 = 0.0
            for k1 in range(0, 6 - i1, 1):
                for k2 in range(0, 7, 1):
                    v10 += x0[i1 + k1 + n10_1,i2 + k2]*mat10[3 + i1 + n10_1,3 + i2,k1,k2]
                
            
            out1[3 + i1 + n10_1,3 + i2] += v10
        
    
    for i1 in range(0, n10_1 + ne10_1, 1):
        for i2 in range(0, ne10_2, 1):
            v10 = 0.0
            for k1 in range(0, 7 - max(0, i1 - n10_1 + 1), 1):
                for k2 in range(0, 6 - i2, 1):
                    v10 += x0[i1 + k1,i2 + k2 + n10_2]*mat10[3 + i1,3 + i2 + n10_2,k1,k2]
                
            
            out1[3 + i1,3 + i2 + n10_2] += v10
        
    
    return