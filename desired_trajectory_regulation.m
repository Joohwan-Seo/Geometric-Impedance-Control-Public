function [Rd, pd, xi_d, dxi_d] = desired_trajectory_regulation(t, g_0)
    w = [0.1,0.1,1]';
    b = 0;
    
    Rd = rodrigues(w,b*t) * g_0(1:3,1:3);
    Rd = [1 0 0;
          0 1 0;
          0 0 1];
    R0 = g_0(1:3,1:3);
%     theta = pi/2 - 0.001;
%     Rd = [1 0 0;
%           0 cos(theta) sin(theta);
%           0 -sin(theta) cos(theta)];

    pd = [-0.5, 0.3, 0.5]';
%     pd = [-0.5, -0.3, 0.2]';
%     pd = g_0(1:3,4);
    
    xi_d = zeros(6,1);
    dxi_d = zeros(6,1);
   
end