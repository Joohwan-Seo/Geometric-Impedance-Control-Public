function [Rd, pd, xi_d, dxi_d] = desired_trajectory3(t,g_0,tf,target)

    if t < 1/5 * tf
        pd = [-0.4, 0.3, 0.2]';
        Rd = [0 -1 0;
              0 0 -1;
              1 0 0];
    elseif t < 2/5 * tf
        pd = [-0.4, -0.3, 0.2]';
        Rd = [1 0 0
              0 0 -1
              0 1 0];
    elseif t < 3/5 * tf
        pd = [-0.6, -0.3, 0.2]';
        Rd = [0 -1 0;
              0 0 -1;
              1 0 0];
    elseif t < 4/5 * tf
        pd = [-0.6, 0.3, 0.2]';
        Rd = [1 0 0
              0 0 -1
              0 1 0];
    else
        pd = [-0.4, 0.3, 0.2]';
        Rd = [0 -1 0;
              0 0 -1;
              1 0 0];
    end
    
    dpd = zeros(3,1);
    ddpd = zeros(3,1);

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