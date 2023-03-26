function [param_x,param_y,param_z,param_w1,param_w2,param_w3] = trajectory_calculator(p_i,p_f,R_i,R_f,ti,tf)    
    A_mat = [[1, ti, ti^2,  ti^3]
             [1, tf, tf^2,  tf^3]
             [0,  1, 2*ti, 3*ti^2]
             [0,  1, 2*tf, 3*tf^2]];
    l_x = [p_i(1), p_f(1), 0, 0]';
    l_y = [p_i(2), p_f(2), 0, 0]';
    l_z = [p_i(3), p_f(3), 0, 0]';
    
    param_x = reshape(inv(A_mat) * l_x, 1, 4);
    param_y = reshape(inv(A_mat) * l_y, 1, 4);
    param_z = reshape(inv(A_mat) * l_z, 1, 4);
    
    init_log = log_map(R_i' * R_i);
    final_log = log_map(R_i' * R_f);

    l_w1 = [init_log(1), final_log(1), 0, 0]';
    l_w2 = [init_log(2), final_log(2), 0, 0]';
    l_w3 = [init_log(3), final_log(3), 0, 0]';
    
    param_w1 = reshape(inv(A_mat) * l_w1, 1, 4);
    param_w2 = reshape(inv(A_mat) * l_w2, 1, 4);
    param_w3 = reshape(inv(A_mat) * l_w3, 1, 4);
end