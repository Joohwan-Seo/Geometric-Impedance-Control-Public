function [Rd, pd, xi_d, dxi_d] = desired_trajectory_regulation(t, g_0)
    w = [0.1,0.1,1]';
    b = 0;
    
    Rd = rodrigues(w,b*t) * g_0(1:3,1:3);
    Rd = [1 0 0;
          0 1 0;
          0 0 1];

    pd = [-0.5, 0.2, 0.3]';
    
    xi_d = zeros(6,1);
    dxi_d = zeros(6,1);
end