function [Rd, pd, xi_d, dxi_d] = desired_trajectory(t, g_0)
    w = [0.1,0.1,1]';
    b = 0;
    
    Rd = rodrigues(w,b*t) * g_0(1:3,1:3);
    
%     Rd = [1 0 0;
%           0 0 1;
%           0 -1 0];
    
%     Rd = g_0(1:3,1:3);

    r1 = 0.20; r2 = 0.20; r3 = 0.05;
    w1 = 1.0; w2 = 1.0; w3 = 0.8;

    pd = [-0.5 - r1 * cos(w1*t); r2 * sin(w2*t); 0.2 + r3 * sin(w3*t)];
    dpd = [r1 * w1 * sin(w1*t); r2 * w2 * cos(w2*t); r3 * w3 * cos(w3*t)];
    ddpd = [r1 * w1^2 * cos(w1*t); -r2 * w2^2 * sin(w2*t); -r3 * w3^2 *sin(w3*t)]; 
    
    wd = Rd' * w * b*t;
    vd = Rd' * dpd;
    
    xi_d = [vd; wd];

    dwd = -hat_map(wd) * Rd' * w + Rd' * w * b;
    dvd = -hat_map(wd) * Rd' * dpd + Rd' * ddpd;

    dxi_d = [dvd; dwd];
end