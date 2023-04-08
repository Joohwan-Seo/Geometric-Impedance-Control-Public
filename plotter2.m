clear; close all; clc;

addpath 'sub_direct'
addpath('results')

obj = 'tracking';
% obj = 'tracking2';
% obj = 'regulation2';
% obj = 'regulation';
%%
if strcmp(obj,'regulation')
    load('results/result_geo2_regulation.mat');
    load('results/result_imp_regulation.mat');
elseif strcmp(obj,'tracking')
    load('results/result_geo2_tracking.mat');
    load('results/result_imp_tracking.mat');
elseif strcmp(obj,'tracking2')
    load('results/result_geo2_tracking2.mat');
    load('results/result_imp_tracking2.mat');
elseif strcmp(obj,'regulation2')
    load('results/result_geo2_regulation2.mat');
    load('results/result_imp_regulation2.mat');
end

%downsample data by factor
down = 10; % default 1

t = result_geo2.t;

N = length(t);
t_ = t(1:end-1);

gd = result_geo2.g_d;

if down > 1
    gd = downsample(gd, down);
    N = ceil(N/down);
    t = downsample(t, down);
    t_ = downsample(t_, down);
    
    result_geo2.S = downsample(result_geo2.S', down)';
    result_geo2.g = downsample(result_geo2.g, down);
    result_geo2.g_d = downsample(result_geo2.g_d, down);
    result_geo2.t = downsample(result_geo2.t, down);
    result_geo2.lyap = downsample(result_geo2.lyap, down);
    result_geo2.potential = downsample(result_geo2.potential, down);
    
    result_imp.S = downsample(result_imp.S', down)';
    result_imp.g = downsample(result_imp.g, down);
    result_imp.g_d = downsample(result_imp.g_d, down);
    result_imp.t = downsample(result_imp.t, down);
    result_imp.lyap = downsample(result_imp.lyap, down);
    result_imp.potential = downsample(result_imp.potential, down);
end

%%
p_geo = zeros(3,N-1);
R_geo = zeros(3,3,N-1);
p_imp = zeros(3,N-1);
R_imp = zeros(3,3,N-1);
p_des = zeros(3,N-1);

for k = 1 : N-1
    p_geo(:,k) = result_geo2.g(k,1:3,4);
    R_geo(:,:,k) = result_geo2.g(k,1:3,1:3);
    p_imp(:,k) = result_imp.g(k,1:3,4);
    R_imp(:,:,k) = result_imp.g(k,1:3,1:3);
    p_des(:,k) = result_geo2.g_d(k,1:3,4);
end

h = 1:N-1/down;

figure(1)
plot3(p_geo(1,h),p_geo(2,h),p_geo(3,h),'r'); hold on; grid on;
plot3(p_imp(1,h),p_imp(2,h),p_imp(3,h),'b--');
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
%%
figure(2)
plot(p_geo(2,h),p_geo(1,h),'r'); hold on; grid on;
plot(p_imp(2,h),p_imp(1,h),'b--');
xlabel('y (m)'); ylabel('x (m)');
nk = 500/down;

for k = 1 : length(h)/nk
    idx = nk*(k-1) + 1;
    scale = 0.05;
    plot([p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,1,idx)],...
         [p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,1,idx)],'r');
    plot([p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,2,idx)],...
         [p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,2,idx)],'b');
    plot([p_geo(2,idx), p_geo(2,idx) + scale * R_geo(2,3,idx)],...
         [p_geo(1,idx), p_geo(1,idx) + scale * R_geo(1,3,idx)],'g');
      
    plot([p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,1,idx)], ...
         [p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,1,idx)],'r');
    plot([p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,2,idx)], ...
         [p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,2,idx)],'b');
    plot([p_imp(2,idx), p_imp(2,idx) + scale * R_imp(2,3,idx)], ...
         [p_imp(1,idx), p_imp(1,idx) + scale * R_imp(1,3,idx)],'g');
end
axis equal
%%
figure(3)
subplot(3,1,1);
plot(t_,gd(1:end-1,1,4),'k:'); hold on; grid on; ylabel('x (m)');
plot(t_,result_geo2.g(1:end-1,1,4),'r'); 
plot(t_,result_imp.g(1:end-1,1,4),'b--');
subplot(3,1,2);
plot(t_,gd(1:end-1,2,4),'k:'); hold on; grid on; ylabel('y (m)');
plot(t_,result_geo2.g(1:end-1,2,4),'r'); 
plot(t_,result_imp.g(1:end-1,2,4),'b--');
subplot(3,1,3);
plot(t_,gd(1:end-1,3,4),'k:'); hold on; grid on; 
plot(t_,result_geo2.g(1:end-1,3,4),'r'); 
plot(t_,result_imp.g(1:end-1,3,4),'b--');
ylabel('z (m)'); xlabel('t (s)');
%%
err_geo = zeros(1,N);
err_imp = zeros(1,N);

err_geo_rot = zeros(1,N);
err_imp_rot = zeros(1,N);

for k = 1 : N
    err_geo(k) = err_fun(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
    err_imp(k) = err_fun(result_imp.g(k,:,:),result_imp.g_d(k,:,:));
    
    err_geo_rot(k) = err_fun_rot(result_geo2.g(k,:,:),result_geo2.g_d(k,:,:));
    err_imp_rot(k) = err_fun_rot(result_imp.g(k,:,:),result_imp.g_d(k,:,:));
end

%%
figure(4)
plot(t_,result_geo2.potential(1:end-1),'r'); hold on; grid on;
plot(t_,result_imp.potential(1:end-1),'b--');
xlabel('t (s)'); ylabel('potential function P')

figure(5)
plot(t_,result_geo2.lyap(1:end-1),'r'); hold on; grid on;
plot(t_,result_imp.lyap(1:end-1),'b--');
xlabel('t (s)'); ylabel('Lyapunov function V');

%%
RMS_p_geo = rms(p_geo - p_des,2);
RMS_p_imp = rms(p_imp - p_des,2);

fprintf('GIC P RMS: %f, %f, %f\n',RMS_p_geo(1),RMS_p_geo(2),RMS_p_geo(3))
fprintf('CIC P RMS: %f, %f, %f\n',RMS_p_imp(1),RMS_p_imp(2),RMS_p_imp(3))

fprintf('GIC ERR RMS: %f, ROT RMS:%f\n', rms(err_geo), rms((err_geo_rot)));
fprintf('CIC ERR RMS: %f, ROT RMS:%f\n', rms(err_imp), rms((err_imp_rot)));

fprintf('GIC Potential RMS: %f\n', rms(result_geo2.potential));
fprintf('CIC Potential RMS: %f\n', rms(result_imp.potential));

fprintf('GIC Lyapunov RMS: %f\n', rms(result_geo2.lyap));
fprintf('CIC Lyapunov RMS: %f\n', rms(result_imp.lyap));

fprintf('\n')
fprintf('===Note===\n')
fprintf('If Lyapunov RMS is ridicuoulsly large, it means that some parts of trajectory are near singular.\n');