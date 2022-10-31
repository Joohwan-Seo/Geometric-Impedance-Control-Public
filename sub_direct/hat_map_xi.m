function mat = hat_map_xi(xi)
    v = xi(1:3); w = xi(4:6);
    v = reshape(v,[3,1]);
    w1 = w(1);
    w2 = w(2);
    w3 = w(3);

    mat_11 = [0, -w3, w2;
              w3, 0, -w1;
              -w2, w1, 0];
    mat = [mat_11, v;
           zeros(1,4)];
end