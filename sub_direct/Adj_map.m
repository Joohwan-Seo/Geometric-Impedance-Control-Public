function mat = Adj_map(input_mat)
    R = input_mat(1:3,1:3);
    p = input_mat(1:3,4);
    
    mat_11 = R;
    mat_12 = hat_map(p) * R;
    mat_21 = zeros(3,3);
    mat_22 = R;
    
    mat = [mat_11, mat_12;
           mat_21, mat_22];
    
end