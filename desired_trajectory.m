function [Rd, pd, xi_d, dxi_d] = desired_trajectory(t, g_0, target)
    w = [0.1,0.1,1]';
    b = 0;
    
    Rd = rodrigues(w,b*t) * g_0(1:3,1:3);

    r1 = 0.15; r2 = 0.15; r3 = 0.10;
    w1 = 2; w2 = 2; w3 = 1;% 
    
    pd = [-0.50 - r1 * cos(w1*t); 0.2 + r2 * sin(w2*t); 0.25 + r3 * sin(w3*t)];
    dpd = [r1 * w1 * sin(w1*t); r2 * w2 * cos(w2*t); r3 * w3 * cos(w3*t)];
    ddpd = [r1 * w1^2 * cos(w1*t); -r2 * w2^2 * sin(w2*t); -r3 * w3^2 *sin(w3*t)];
    
    Rd_dot = zeros(3,3);
    Rd_ddot = zeros(3,3);
    if strcmp(target, 'geo')
        wd = vee_map(Rd' * Rd_dot);
        vd = Rd' * dpd;
        dwd = vee_map(Rd_dot' * Rd_dot + Rd' * Rd_ddot);
        dvd = Rd_dot' * dpd + Rd' * ddpd;
    elseif strcmp(target, 'imp')
        wd = vee_map(Rd_dot * Rd');
        vd = dpd;
        dwd = vee_map(Rd_dot * Rd_dot' + Rd_ddot * Rd');
        dvd = ddpd;
    end        
    
    xi_d = [vd; wd];
    dxi_d = [dvd; dwd];
end