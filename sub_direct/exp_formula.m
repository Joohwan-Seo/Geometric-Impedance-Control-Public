function mat = exp_formula(xi, q)
    v = xi(1:3); w = xi(4:6);
    
    w_hat = hat_map(w);
    exp_w = rodrigues(w,q);
   
    mat_11 = exp_w;
    mat_12 = (eye(3) - exp_w) * w_hat * v + w * transpose(w) * v * q;
    
    mat = [mat_11, mat_12;
           zeros(1,3), 1];
end