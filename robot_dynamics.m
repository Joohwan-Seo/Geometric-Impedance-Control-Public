function dS = robot_dynamics(t,q,T)
    q1 = q(1); q2 = q(2); q3 = q(3);
    q4 = q(4); q5 = q(5); q6 = q(6);
    dq1 = q(7); dq2 = q(8); dq3 = q(9);
    dq4 = q(10); dq5 = q(11); dq6 = q(12);
    
    dS = zeros(12,1);
    dq = [dq1, dq2, dq3, dq4, dq5, dq6]';
    M = M_fun(q1,q2,q3,q4,q5,q6);
    C = C_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    G = G_fun(q1,q2,q3,q4,q5,q6);
%%
    dS(1) = dq1;
    dS(2) = dq2;
    dS(3) = dq3;
    dS(4) = dq4;
    dS(5) = dq5;
    dS(6) = dq6;
    
    J = Je_fun(q1,q2,q3,q4,q5,q6);


    Fx = 0; Fy = 0; Fz = 0;
    tx = 0; ty = 0; tz = 0;
    
    Fe = [Fx; Fy; Fz; tx; ty; tz];
    
    Te = J' * Fe;

    ddq = inv(M) * (-C*dq - G' + T + Te);
    
    dS(7) = ddq(1);
    dS(8) = ddq(2);
    dS(9) = ddq(3);
    dS(10) = ddq(4);
    dS(11) = ddq(5);
    dS(12) = ddq(6);
end