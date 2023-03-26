clear; close all; clc;
%%
% Comparison for two geometric impedance controllers
%%
addpath 'sub_direct'
addpath 'toTIKZ'

obj = 'tracking'; % or regulation
% obj = 'tracking2';
% obj = 'regulation2';
% obj = 'regulation';
%%
export_to_tikz = false;
%%
if strcmp(obj,'regulation')
    load('result_geo_regulation.mat');
    load('result_geo2_regulation.mat');
elseif strcmp(obj,'tracking')
    load('result_geo_tracking.mat');
    load('result_geo2_tracking.mat');
elseif strcmp(obj,'tracking2')
    load('result_geo_tracking2.mat');
    load('result_geo2_tracking2.mat');
elseif strcmp(obj,'regulation2')
    load('result_geo_regulation2.mat');
    load('result_geo2_regulation2.mat');
end

%downsample data by factor
down = 1; % default 1

t = result_geo.t;

N = length(t);
t_ = t(1:end-1);

gd = result_geo.g_d;

if down > 1
    gd = downsample(gd, down);
    N = ceil(N/down);
    t = downsample(t, down);
    t_ = downsample(t_, down);
    result_geo.S = downsample(result_geo.S', down)';
    result_geo.g = downsample(result_geo.g, down);
    result_geo.g_d = downsample(result_geo.g_d, down);
    result_geo.t = downsample(result_geo.t, down);
    result_geo2.S = downsample(result_geo2.S', down)';
    result_geo2.g = downsample(result_geo2.g, down);
    result_geo2.g_d = downsample(result_geo2.g_d, down);
    result_geo2.t = downsample(result_geo2.t, down);
    
    result_geo.lyap = downsample(result_geo.lyap, down);
    result_geo2.lyap = downsample(result_geo2.lyap, down);
end

%%
p_geo = zeros(3,N-1);
R_geo = zeros(3,3,N-1);
p_imp = zeros(3,N-1);
R_imp = zeros(3,3,N-1);
p_des = zeros(3,N-1);

for k = 1 : N-1
    p_geo(:,k) = result_geo.g(k,1:3,4);
    R_geo(:,:,k) = result_geo.g(k,1:3,1:3);
    p_imp(:,k) = result_geo2.g(k,1:3,4);
    R_imp(:,:,k) = result_geo2.g(k,1:3,1:3);
    p_des(:,k) = result_geo.g_d(k,1:3,4);
end

h = 1:N-1/down;

figure(1)
plot3(p_geo(1,h),p_geo(2,h),p_geo(3,h),'k'); hold on; grid on;
plot3(p_imp(1,h),p_imp(2,h),p_imp(3,h),'m');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
nk = 500/down;

for k = 1 : length(h)/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,1,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,1,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,1,idx)],'r');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,2,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,2,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,2,idx)],'b');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,3,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,3,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,3,idx)],'g');
      
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,1,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,1,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(3,1,idx)],'r');
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,2,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,2,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(3,2,idx)],'b');
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,3,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,3,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(3,3,idx)],'g');
end
axis equal

%export fig as tex-file
if export_to_tikz
    warning('off')
    name = ['pics\trajectory_3d_', obj , '.tex'];
    matlab2tikz(name)
    warning('on')
end
%%
figure(2)
subplot(3,1,1);
plot(t_,gd(1:end-1,1,4),'k:'); hold on; grid on; ylabel('x (m)');
plot(t_,result_geo.g(1:end-1,1,4),'r'); 
plot(t_,result_geo2.g(1:end-1,1,4),'b--');
subplot(3,1,2);
plot(t_,gd(1:end-1,2,4),'k:'); hold on; grid on; ylabel('y (m)');
plot(t_,result_geo.g(1:end-1,2,4),'r'); 
plot(t_,result_geo2.g(1:end-1,2,4),'b--');
subplot(3,1,3);
plot(t_,gd(1:end-1,3,4),'k:'); hold on; grid on; 
plot(t_,result_geo.g(1:end-1,3,4),'r'); 
plot(t_,result_geo2.g(1:end-1,3,4),'b--');
ylabel('z (m)'); xlabel('t (s)');

%export fig as tex-file
if export_to_tikz
    warning('off')
    name = ['pics\comparison_xyz_', obj , '.tex'];
    matlab2tikz(name)
    warning('on')
end
%%
err_geo = zeros(1,N);
err_imp = zeros(1,N);

err_geo_rot = zeros(1,N);
err_imp_rot = zeros(1,N);

for k = 1 : N
    err_geo(k) = err_fun(result_geo.g(k,:,:),result_geo.g_d(k,:,:));
    err_imp(k) = err_fun(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
    
    err_geo_rot(k) = err_fun_rot(result_geo.g(k,:,:),result_geo.g_d(k,:,:));
    err_imp_rot(k) = err_fun_rot(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
end

%%
figure(3)
plot(t_,err_geo(1:end-1),'r'); hold on; grid on;
plot(t_,err_imp(1:end-1),'b--');
xlabel('t (s)'); ylabel('Error function \Psi')

%export fig as tex-file
if export_to_tikz
    warning('off')
    name = ['pics\error_func_', obj , '.tex']
    matlab2tikz(name)
    warning('on')
end

figure(4)
plot(t_,result_geo.lyap(1:end-1),'r'); hold on; grid on;
plot(t_,result_geo2.lyap(1:end-1),'b--');
xlabel('t (s)'); ylabel('Dynamic Cost \Phi');

if export_to_tikz
    warning('off')
    name = ['pics\dynamic_cost_', obj , '.tex'];
    matlab2tikz(name)
    warning('on')
end
% ylim([0,1])
%%
RMS_p_geo = rms(p_geo - p_des,2);
RMS_p_imp = rms(p_imp - p_des,2);

fprintf('GIC P RMS: %f, %f, %f\n',RMS_p_geo(1),RMS_p_geo(2),RMS_p_geo(3))
fprintf('GIC2 P RMS: %f, %f, %f\n',RMS_p_imp(1),RMS_p_imp(2),RMS_p_imp(3))

fprintf('GIC ERR RMS: %f, ROT RMS:%f\n', rms(err_geo), rms(err_geo_rot));
fprintf('GIC2 ERR RMS: %f, ROT RMS:%f\n', rms(err_imp), rms(err_imp_rot));

fprintf('GIC Dynamic cost RMS: %f\n', rms(result_geo.lyap));
fprintf('GIC2 Dynamic cost RMS: %f\n', rms(result_geo2.lyap));

%% Testing for rightful rotation matrix
det_imp = zeros(1,N-1);
det_geo = zeros(1,N-1);

for k = 1 : N-1
    det_imp(k) = det(R_imp(:,:,k));
    det_geo(k) = det(R_geo(:,:,k));
end

figure(3)
plot(t_,err_geo(1:end-1),'r'); hold on; grid on;
plot(t_,err_imp(1:end-1),'b--');
xlabel('t (s)'); ylabel('Error function \Psi')


