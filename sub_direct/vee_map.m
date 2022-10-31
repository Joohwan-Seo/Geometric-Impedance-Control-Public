function vec = vee_map(mat)
    vec = zeros(3,1);
    
    vec(3) = -mat(1,2);
    vec(1) = -mat(2,3);
    vec(2) = mat(1,3);  
    
end