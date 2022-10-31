clear; close all; clc;

addpath 'sub_direct'
addpath 'toTIKZ'
%%
export_to_tikz = false;
%%
% load('result_geo_regulation.mat');
% load('result_imp_regulation.mat');

load('result_geo_tracking.mat');
load('result_imp_tracking.mat');

%downsample data by factor
down = 10; % default 1

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
    result_imp.S = downsample(result_imp.S', down)';
    result_imp.g = downsample(result_imp.g, down);
    result_imp.g_d = downsample(result_imp.g_d, down);
    result_imp.t = downsample(result_imp.t, down);
end


figure(1)
subplot(3,1,1);
plot(t_,gd(1:end-1,1,4),'k:'); hold on; grid on; ylabel('x (m)');
plot(t_,result_geo.g(1:end-1,1,4),'r'); 
plot(t_,result_imp.g(1:end-1,1,4),'b--');
subplot(3,1,2);
plot(t_,gd(1:end-1,2,4),'k:'); hold on; grid on; ylabel('y (m)');
plot(t_,result_geo.g(1:end-1,2,4),'r'); 
plot(t_,result_imp.g(1:end-1,2,4),'b--');
subplot(3,1,3);
plot(t_,gd(1:end-1,3,4),'k:'); hold on; grid on; ylabel('z (m)'); xlabel('t (s)');
plot(t_,result_geo.g(1:end-1,3,4),'r'); 
plot(t_,result_imp.g(1:end-1,3,4),'b--');

%export fig as tex-file
if export_to_tikz
    warning('off')
    matlab2tikz('pics\comparison_xyz.tex')
    warning('on')
end


%%
p_geo = zeros(3,N-1);
R_geo = zeros(3,3,N-1);
p_imp = zeros(3,N-1);
R_imp = zeros(3,3,N-1);
for k = 1 : N-1
    p_geo(:,k) = result_geo.g(k,1:3,4);
    R_geo(:,:,k) = result_geo.g(k,1:3,1:3);
    p_imp(:,k) = result_imp.g(k,1:3,4);
    R_imp(:,:,k) = result_imp.g(k,1:3,1:3);
end



h = 1:10000/down;
figure(2)
plot3(p_geo(1,h),p_geo(2,h),p_geo(3,h),'k'); hold on; grid on;
plot3(p_imp(1,h),p_imp(2,h),p_imp(3,h),'m');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
nk = 500/down;

for k = 1 : length(h)/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,1,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(1,2,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(1,3,idx)],'r');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(2,1,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,2,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(2,3,idx)],'b');
    plot3([p_geo(1,idx), p_geo(1,idx) + scale * R_geo(3,1,idx)],...
          [p_geo(2,idx), p_geo(2,idx) + scale * R_geo(3,2,idx)],...
          [p_geo(3,idx), p_geo(3,idx) + scale * R_geo(3,3,idx)],'g');
      
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,1,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(1,2,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(1,3,idx)],'r');
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(2,1,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,2,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(2,3,idx)],'b');
    plot3([p_imp(1,idx), p_imp(1,idx) + scale * R_imp(3,1,idx)], ...
          [p_imp(2,idx), p_imp(2,idx) + scale * R_imp(3,2,idx)], ...
          [p_imp(3,idx), p_imp(3,idx) + scale * R_imp(3,3,idx)],'g');
end

%%
err_geo = zeros(1,N);
err_imp = zeros(1,N);
for k = 1 : N
    err_geo(k) = err_fun(result_geo.g(k,:,:),result_geo.g_d(k,:,:));
    err_imp(k) = err_fun(result_imp.g(k,:,:),result_imp.g_d(k,:,:));
end

%export fig as tex-file
if export_to_tikz
    warning('off')
    matlab2tikz('pics\trajectory_3d.tex')
    warning('on')
end

figure(3)
plot(t_,err_geo(1:end-1),'r'); hold on; grid on;
plot(t_,err_imp(1:end-1),'b--');
xlabel('t (s)'); ylabel('Error function \Psi')


%% Testing for rightful rotation matrix
det_imp = zeros(1,N-1);
det_geo = zeros(1,N-1);
for k = 1 : N-1
    det_imp(k) = det(R_imp(:,:,k));
    det_geo(k) = det(R_geo(:,:,k));
end

%export fig as tex-file
if export_to_tikz
    warning('off')
    matlab2tikz('pics\error_func.tex')
    warning('on')
end