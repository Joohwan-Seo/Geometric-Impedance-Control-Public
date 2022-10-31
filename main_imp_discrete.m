clear; close all; clc;
%%
saving = false;

addpath('sub_direct');
obj = 'tracking'; % or regulation
% obj = 'regulation'; % or regulation
%% initialization
if strcmp(obj,'tracking')
    t_end = 20;
elseif strcmp(obj,'regulation')
    t_end = 10;
end

t = 0 : 0.001 : t_end;
N = length(t);

S = zeros(12,N);
% S(1:6,1) = [0.1, -0.3, 0.4, 0.2, 0.5, 0.2];

g_0 = g_st_fun(0,0,0,0,0,0);

T_arr = zeros(6,N);
T1_arr = zeros(6,N);
g_se_arr = zeros(N,4,4);
pd_arr = zeros(3,N);
g_d_arr = zeros(N,4,4);
F_arr = zeros(6,N);

%% simulation
for k = 1 : N-1
    %% input calculation
    gain_t = 100;
    gain_o = 100;
    kt1 = gain_t; kt2 = gain_t; kt3 = gain_t;
    ko1 = gain_o; ko2 = gain_o; ko3 = gain_o;
    
    kv = 50;
    
    Kt = diag([kt1,kt2,kt3]);
    Ko = diag([ko1,ko2,ko3]);

    T_max = [150, 100, 100, 100, 100, 100]';
    
    Kg = blkdiag(Kt,Ko);
    Kxi = kv * eye(6);
    
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
    Je = Je_fun(q1,q2,q3,q4,q5,q6);
    Je_dot = Je_dot_fun(q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6);
    G = G_fun(q1,q2,q3,q4,q5,q6);
    g_se = g_st_fun(q1,q2,q3,q4,q5,q6);
    
    g_se_arr(k,:,:) = g_se;
    
    R = g_se(1:3,1:3); p = g_se(1:3,4);
    xi = Je * dq;
    v = xi(1:3); w = xi(4:6);
    
    %%
    if strcmp(obj,'tracking')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory_spatial(t(k),g_0);
    elseif strcmp(obj,'regulation')
        [Rd,pd,xi_d,dxi_d] = desired_trajectory_regulation(t(k),g_0);
    end
    
    pd_arr(:,k) = pd;
    g_d_arr(k,:,:) = [Rd, pd;
                     zeros(1,3),1];
    
    vd = xi_d(1:3); wd = xi_d(4:6);
    g_ed = [R'*Rd, R'*(p - pd);
            zeros(1,3),1];
    R_ed = g_ed(1:3,1:3);
    p_ed = g_ed(1:3,4);

    %%
    rc1 = R(:,1); rc2 = R(:,2); rc3 = R(:,3);
    rd1 = Rd(:,1); rd2 = Rd(:,2); rd3 = Rd(:,3);
%     rot_err = 0.5 * (cross(rc1,rd1) + cross(rc2, rd2) + cross(rc3, rd3));
    rot_err = 0.5* vee_map(Rd'*R - R'*Rd);

    %%
    R_ed_dot = -hat_map(w)*R'*Rd + R'*Rd*hat_map(wd);
    p_ed_dot = -hat_map(w)*R'*(p - pd) + v - R'*Rd*vd;
    Ad_ged_dot = [R_ed_dot, hat_map(p_ed_dot)*R_ed + hat_map(p_ed)*R_ed_dot;
                  zeros(3,3), R_ed_dot];
    
%     P = Ad_ged_dot * xi_d + Adj_map(g_ed) * dxi_d;
%     
%     xi_d_t = Adj_map(g_ed) * xi_d;

    e_x = [p - pd; rot_err];
    wd_t = R_ed * wd;
    xi_d = [vd ; wd_t];
    de_x = xi - xi_d;
    ddwd = R'*-(hat_map(w))*Rd*wd;
    dde_x = [dxi_d(1:3,1);ddwd];
    %%
    Je_inv = inv(M) * Je' * pinv(Je * inv(M) * Je'); %correct this with pseudo inverse
    M_tilde = Je_inv' * M * Je_inv;
    C_tilde = (Je_inv)' * (C - M * Je_inv * Je_dot) * Je_inv; 
    
    Ko = Rd * Ko * Rd';
    Kt = Rd * Kt * Rd';
    Kxi = blkdiag(Rd,Rd) * Kxi * blkdiag(Rd', Rd');
    Kg = blkdiag(Kt,Ko);
    
    if abs(det(Je)) < 0.01
        F = - Kg * e_x - Kxi * de_x;
    else
        F = C_tilde * xi - Kg * e_x - Kxi * de_x + M_tilde * dxi_d;
    end
    
    T1 = Je'* F;
    T = T1 + G';
    T_arr(:,k) = T;
    T1_arr(:,k) = T1;
    F_arr(:,k) = F;
    
    initVal = S(:,k);
    [Time,S_] = ode15s(@(Time,S_) robot_dynamics(Time,S_,T), [t(k),t(k+1)], initVal ,' ');
    [Nn,~] = size(S_);
    S(:,k+1) = S_(Nn,:);
    
    if mod(k,1000) == 0
        disp(k/1000)
    end
end

%%
t_ = t(1:end-1);

figure
subplot(3,1,1);
plot(t_,g_se_arr(1:end-1,1,4));
subplot(3,1,2);
plot(t_,g_se_arr(1:end-1,2,4));
subplot(3,1,3);
plot(t_,g_se_arr(1:end-1,3,4));

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
    plot3([p(1,idx), p(1,idx) + scale * R(1,1,idx)],[p(2,idx), p(2,idx) + scale * R(1,2,idx)],...
          [p(3,idx), p(3,idx) + scale * R(1,3,idx)],'r');
    plot3([p(1,idx), p(1,idx) + scale * R(2,1,idx)],[p(2,idx), p(2,idx) + scale * R(2,2,idx)],...
          [p(3,idx), p(3,idx) + scale * R(2,3,idx)],'b');
    plot3([p(1,idx), p(1,idx) + scale * R(3,1,idx)],[p(2,idx), p(2,idx) + scale * R(3,2,idx)],...
          [p(3,idx), p(3,idx) + scale * R(3,3,idx)],'g');
end
%%
if saving == true
result_imp.S = S;
result_imp.g = g_se_arr;
result_imp.g_d = g_d_arr;
result_imp.T = T;
result_imp.T1 = T1;
result_imp.t = t;

    if strcmp(obj,"tracking")
        save('result_imp_tracking.mat','result_imp');
    elseif strcmp(obj,'regulation')
        save('result_imp_regulation.mat','result_imp');
    end
end

%%
T_rms = rms(T,2);
T1_rms = rms(T1,2);
F_rms = rms(F,2);
pos_rms = rms(p - pd_arr,2);

disp('Impedance T rms')
disp(T_rms')
disp(F_rms');

disp('Impedance pos_rms');
disp(pos_rms)