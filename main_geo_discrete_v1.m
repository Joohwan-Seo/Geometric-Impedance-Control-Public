clear; close all; clc;
%%
% This control law is based on the 'geometric impedance control law version 1'
% Equation (32)
%%
saving = true;

addpath('sub_direct')
addpath('results')

obj = 'tracking'; 
% obj = 'tracking2';
% obj = 'regulation2';
% obj = 'regulation';
%% initialization
if strcmp(obj,'tracking')
    t_end = 10;
elseif strcmp(obj,'regulation')
    t_end = 5;
elseif strcmp(obj,'tracking2')
    t_end = 5;
elseif strcmp(obj,'regulation2')
    t_end = 15;
end
t = 0 : 0.001 : t_end;
N = length(t);

S = zeros(12,N);
if strcmp(obj,'tracking')
   S(1:6,1) = [0.2, -0.5, 0.4, 0.6, -0.5, 0.2];
elseif strcmp(obj,'tracking2')
   S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, 0];

   g_0 = g_st_fun(S(1,1),S(2,1),S(3,1),S(4,1),S(5,1),S(6,1));
   p_i = [-0.5, -0.3, 0.2]';
   p_f = [-0.5, 0.3, 0.2]';
   R_i = g_0(1:3,1:3);
   R_f = [0 -1 0;
          0 0 -1;
          1 0 0];
   [param_x,param_y,param_z,param_w1,param_w2,param_w3] = trajectory_calculator(p_i,p_f,R_i,R_f,t(1),t(end));
   param_p = [param_x;  param_y;  param_z];
   param_w = [param_w1; param_w2; param_w3];

   S(1:6,1) = [0.4, -pi/6, 0.4, 0.6, -0.5, 0.2];
elseif strcmp(obj,'regulation')
    S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, 0];
elseif strcmp(obj,'regulation2')
    S(1:6,1) = [0.1721, -1.0447, 1.6729, -0.6282, 0.1721, 0];
end

g_0 = g_st_fun(0,0,0,0,0,0);

T_arr = zeros(6,N);
T1_arr = zeros(6,N);
g_se_arr = zeros(N,4,4);
pd_arr = zeros(3,N);
g_d_arr = zeros(N,4,4);
F_arr = zeros(6,N);
sing_arr = zeros(1,N);
err_arr = zeros(1,N);
err_rot_arr = zeros(1,N);

lyap_arr = zeros(1,N);
potential_arr = zeros(1,N);

%% simulation
for k = 1 : N-1
    %% input calculation  
    if strcmp(obj,'tracking')
        kp1 = 200; kp2 = 60; kp3 = 80;
        kr1 = 10; kr2 = 30; kr3 = 100;
    else
        kp1 = 100; kp2 = 100; kp3 = 100;
        kr1 = 100; kr2 = 100; kr3 = 100;
    end
    
    kd = 50;
    
    Kp = diag([kp1,kp2,kp3]);
    KR = diag([kr1,kr2,kr3]);
    
    Kg = blkdiag(Kp,KR);
    Kd = kd * eye(6);    
    %%
    q1 = S(1,k); q2 = S(2,k); q3 = S(3,k); 
    q4 = S(4,k); q5 = S(5,k); q6 = S(6,k);
    dq1 = S(7,k); dq2 = S(8,k); dq3 = S(9,k); 
    dq4 = S(10,k); dq5 = S(11,k); dq6 = S(12,k);
    
    q = [q1, q2, q3, q4, q5, q6]';
    dq = [dq1, dq2, dq3, dq4, dq5, dq6]';
    %%
    M = M_fun(q1,q2,q3,q4,q5,q6);
    C = C_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    Jb = Jb_fun(q1,q2,q3,q4,q5,q6);
    Jb_dot = Jb_dot_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    G = G_fun(q1,q2,q3,q4,q5,q6);
    g_se = g_st_fun(q1,q2,q3,q4,q5,q6);
    
    g_se_arr(k,:,:) = g_se;
    
    R = g_se(1:3,1:3); p = g_se(1:3,4);
    V_b = Jb * dq; %Body velocity
    v = V_b(1:3); w = V_b(4:6);
    %%
    if strcmp(obj,'tracking')
        [Rd,pd,V_d,dV_d] = desired_trajectory(t(k),g_0,'geo');
    elseif strcmp(obj,'regulation')
        [Rd,pd,V_d,dV_d] = desired_trajectory_regulation(t(k),g_0);
    elseif strcmp(obj,'tracking2')
        [Rd,pd,V_d,dV_d] = desired_trajectory2(t(k),g_0,param_p, param_w,'geo');
    elseif strcmp(obj,'regulation2')
        [Rd,pd,V_d,dV_d] = desired_trajectory_regulation2(t(k),g_0,t(end),'geo');
    end
    
    vd = V_d(1:3); wd = V_d(4:6);
    g_ed = [R'*Rd, R'*(p - pd);
            zeros(1,3),1];
        
    R_ed = g_ed(1:3,1:3);
    p_ed = g_ed(1:3,4);
    
    pd_arr(:,k) = pd;
    g_d_arr(k,:,:) = [Rd, pd;
                      zeros(1,3),1];

    f_R = vee_map(KR *Rd'*R - R'*Rd * KR); 
    f_p = R'* (Rd * Kp * Rd') * (p - pd);

    V_d_star = Adj_map(g_ed) * V_d;
    
    f_g = [f_p; f_R];
    e_V = V_b - V_d_star;
    %%
    R_ed_dot = -hat_map(w)*R'*Rd + R'*Rd*hat_map(wd);
    p_ed_dot = -hat_map(w)*R'*(p - pd) + v - R'*Rd*vd;
    Ad_ged_dot = [R_ed_dot, hat_map(p_ed_dot)*R_ed + hat_map(p_ed)*R_ed_dot;
                  zeros(3,3), R_ed_dot];
    
    dV_d_star = Ad_ged_dot * V_d + Adj_map(g_ed) * dV_d;
    %%
    Jb_inv = inv(M) * Jb' * pinv(Jb * inv(M) * Jb'); 
    M_tilde = Jb_inv' * M * Jb_inv;
    C_tilde = (Jb_inv)' * (C - M * Jb_inv * Jb_dot) * Jb_inv; 
    
    if abs(det(Jb)) < 0.01
        sing_arr(k) = 1;
        F = - f_g - Kd * e_V;
    else
        F = C_tilde * V_d_star -  f_g - Kd * e_V + M_tilde * dV_d_star;
    end
    
    
    T1 = Jb'* F;
    T = T1 + G';
    T_arr(:,k) = T;
    T1_arr(:,k) = T1;
    F_arr(:,k) = F;
    
    initVal = S(:,k);
    [Time,S_] = ode15s(@(Time,S_) robot_dynamics(Time,S_,T), [t(k),t(k+1)], initVal ,' ');
    [Nn,~] = size(S_);
    S(:,k+1) = S_(Nn,:);
 
    potential = potential_fun(g_se, g_d_arr(k,:,:),Kp,KR);
    kinetic = 1/2 * e_V' * M_tilde * e_V;
    
    err_arr(k) = err_fun(g_se_arr(k,:,:),g_d_arr(k,:,:));
    err_rot_arr(k) = err_fun_rot(g_se_arr(k,:,:),g_d_arr(k,:,:)); 

    lyap_arr(k) = potential + kinetic;
    potential_arr(k) = potential;

    if mod(k,1000) == 0
        disp(k/1000)
    end
end

%%
t_ = t(1:end-1);

figure
subplot(3,1,1);
plot(t_,g_se_arr(1:end-1,1,4)); hold on; plot(t_,pd_arr(1,1:end-1))
subplot(3,1,2);
plot(t_,g_se_arr(1:end-1,2,4)); hold on; plot(t_,pd_arr(2,1:end-1))
subplot(3,1,3);
plot(t_,g_se_arr(1:end-1,3,4)); hold on; plot(t_,pd_arr(3,1:end-1))

%%
p = zeros(3,N);
R = zeros(3,3,N);
for k = 1 : N
    p(:,k) = g_se_arr(k,1:3,4);
    R(:,:,k) = g_se_arr(k,1:3,1:3);
end

figure
plot3(g_se_arr(1:end-1,1,4),g_se_arr(1:end-1,2,4),g_se_arr(1:end-1,3,4)); hold on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

nk = 500;
for k = 1 : N/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot3([p(1,idx), p(1,idx) + scale * R(1,1,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,1,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,1,idx)],'r');
    plot3([p(1,idx), p(1,idx) + scale * R(1,2,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,2,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,2,idx)],'b');
    plot3([p(1,idx), p(1,idx) + scale * R(1,3,idx)],...
          [p(2,idx), p(2,idx) + scale * R(2,3,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,3,idx)],'g');
end
axis equal
%%
if saving == true
result_geo1.S = S;
result_geo1.g = g_se_arr;
result_geo1.g_d = g_d_arr;
result_geo1.T = T;
result_geo1.T1 = T1;
result_geo1.t = t;
result_geo1.potential = potential_arr;
result_geo1.lyap = lyap_arr;

    if strcmp(obj,"tracking")
        save('results/result_geo1_tracking.mat','result_geo1');
    elseif strcmp(obj,'regulation')
        save('results/result_geo1_regulation.mat','result_geo1');
    elseif strcmp(obj,'tracking2')
        save('results/result_geo1_tracking2.mat','result_geo1');
    elseif strcmp(obj,'regulation2')
        save('results/result_geo1_regulation2.mat','result_geo1');
    end
end