function [Rd, pd, xi_d, dxi_d] = desired_trajectory2(t,g_0,param_p, param_w,target)

    pd = param_p(:,1) + param_p(:,2) * t + param_p(:,3) * t^2 + param_p(:,4) * t^3;
    wt = param_w(:,1) + param_w(:,2) * t + param_w(:,3) * t^2 + param_w(:,4) * t^3;
    
    dpd = param_p(:,2) + 2 * param_p(:,3) * t + 3 * param_p(:,4) * t^2;
    dwt = param_w(:,2) + 2 * param_w(:,3) * t + 3 * param_w(:,4) * t^2;
    
    ddpd = 2 * param_p(:,3) + 6 * param_p(:,4) * t;
    ddwt = 2 * param_w(:,3) + 6 * param_w(:,4) * t; 
    
    R0 = [1 0 0;
          0 0 -1;
          0 1 0];
%     R0 = [1 0 0;
%           0 1 0;
%           0 0 1];
%     R0 = g_0(1:3,1:3);
    
    Rt = expm(hat_map(wt));
    Rd = R0 * Rt;
    Rd_dot = R0 * hat_map(dwt) * Rt;
    Rd_ddot = R0 * (hat_map(ddwt)*Rt + hat_map(dwt) * hat_map(dwt) * Rt);
    
    
%     Rd_dot = zeros(3,3);
%     Rd_ddot = zeros(3,3);
    
    if strcmp(target, 'geo')
        wd = vee_map(Rd' * Rd_dot);
        vd = Rd' * dpd;
        dwd = vee_map(Rd_dot' * Rd_dot + Rd' * Rd_ddot);
        dvd = Rd_dot' * dpd + Rd' * ddpd;
    elseif strcmp(target, 'imp')
        wd = vee_map(Rd_dot * Rd');
%         vd = dpd - Rd_dot * Rd' * pd;
        vd = dpd;
        dwd = vee_map(Rd_dot * Rd_dot' + Rd_ddot * Rd');
%         dvd = ddpd - Rd_ddot * Rd' * pd - Rd_dot * Rd_dot' * pd - Rd_dot * Rd' * dpd;
        dvd = ddpd;
    end        
    
    xi_d = [vd; wd];
    dxi_d = [dvd; dwd];
end