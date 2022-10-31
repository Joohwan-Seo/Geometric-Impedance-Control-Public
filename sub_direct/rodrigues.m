function rot_mat = rodrigues(axis,theta)
    w = axis;
    w_abs = norm(w,2);

    w_hat = hat_map(w);
    n = length(axis);
    
    rot_mat = eye(n) + sin(w_abs * theta)/w_abs * w_hat + (1 - cos(w_abs * theta))/w_abs^2 * w_hat^2;

end