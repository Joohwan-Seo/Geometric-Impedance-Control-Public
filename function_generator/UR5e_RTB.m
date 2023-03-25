clear; close all; clc;
target = 'matlab';
%%
if strcmp(target,'xml')
alpha = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.163, 0, 0, 0.134, 0.1, 0.0771]; % d is given in z-direction
a = [0, -0.425, -0.392, 0, 0, 0]; % a is given in x-direction
m = [3.7, 8.393, 2.275, 1.219, 1.219, 0.1889];

p_c1 = [0, 0, 0];
p_c2 = [0.2125, 0, 0.138];
p_c3 = [0.196, 0, 0.007];
p_c4 = [0, 0, 0];
p_c5 = [0, 0, 0];
p_c6 = [0, 0, 0];

In_1 = diag([0.010267, 0.010267, 0.00660]);
In_2 = diag([0.015107, 0.13389, 0.13389]);
In_3 = diag([0.004095, 0.031178, 0.031178]);
In_4 = diag([0.0025599, 0.0021942, 0.0025599]);
In_5 = diag([0.0025599, 0.0025599, 0.0021942]);
In_6 = diag([0.00009804, 0.00009804, 0.00013210]);

end
%% From robosuite XML files
% In_1 = diag([0.010267, 0.010267, 0.00660]);
% In_2 = diag([0.015107, 0.13389, 0.13389]);
% In_3 = diag([0.004095, 0.031178, 0.031178]);
% In_4 = diag([0.0025599, 0.0021942, 0.0025599]);
% In_5 = diag([0.0025599, 0.0025599, 0.0021942]);
% In_6 = diag([0.00009804, 0.00009804, 0.00013210]);
%% From matlab UR5 moment of inertia file
if strcmp(target,'matlab')
alpha = [pi/2 0 0 pi/2 -pi/2 0];
d = [0.163, 0, 0, 0.134, 0.1, 0.0771]; % d is given in z-direction
a = [0, -0.425, -0.39225, 0, 0, 0]; % a is given in x-direction
m = [3.7, 8.393, 2.33, 1.219, 1.219, 0.1889];

p_c1 = [0, -0.02561, 0.00193];
p_c2 = [0.2125, 0, 0.11336];
p_c3 = [0.196, 0, 0.0265];
p_c4 = [0, -0.0018, 0.01634];
p_c5 = [0, -0.0018, 0.01634];
p_c6 = [0, 0, -0.001159];

In_1 = [0.010267, 0.00660, 0.010267];
In_2 = [0.0151, 0.8849, 0.8849];
In_3 = [0.004095, 0.1916, 0.1916];
In_4 = [0.1112, 0.2194, 0.1112];
In_5 = [0.1112, 0.2194, 0.1112];
In_6 = [0.0171, 0.0171, 0.0338];
end
%%

L(1) = Revolute('d', d(1), ...   
    'a', a(1), ...               
    'alpha', alpha(1), ...        
    'I', [In_1(1), In_1(2), In_1(3), 0, 0, 0], ... 
    'r', p_c1, ...       
    'm', m(1)); 

L(2) = Revolute('d', d(2), ...   
    'a', a(2), ...               
    'alpha', alpha(2), ...        
    'I', [In_2(1), In_2(2), In_2(3), 0, 0, 0], ... 
    'r', p_c2, ...       
    'm', m(2)); 

L(3) = Revolute('d', d(3), ...   
    'a', a(3), ...               
    'alpha', alpha(3), ...        
    'I', [In_3(1), In_3(2), In_3(3), 0, 0, 0], ... 
    'r', p_c1, ...       
    'm', m(1)); 

L(4) = Revolute('d', d(4), ...   
    'a', a(4), ...               
    'alpha', alpha(4), ...        
    'I', [In_4(1), In_4(2), In_4(3), 0, 0, 0], ... 
    'r', p_c4, ...       
    'm', m(4)); 

L(5) = Revolute('d', d(5), ...   
    'a', a(5), ...               
    'alpha', alpha(5), ...        
    'I', [In_5(1), In_5(2), In_5(3), 0, 0, 0], ... 
    'r', p_c5, ...       
    'm', m(5)); 

L(6) = Revolute('d', d(6), ...   
    'a', a(6), ...               
    'alpha', alpha(6), ...        
    'I', [In_6(1), In_6(2), In_6(3), 0, 0, 0], ... 
    'r', p_c6, ...       
    'm', m(6)); 

UR5e = SerialLink(L, 'name', 'UR5e');

UR5e.sym()
%%
syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real

q = [q1, q2, q3, q4, q5, q6];
dq = [dq1, dq2, dq3,dq4,dq5,dq6];

tic
% M = UR5e.inertia(q);
Jb = UR5e.jacobe(q);
% C = UR5e.coriolis(q,dq);
% G = UR5e.gravload(q);
g_st_ = UR5e.fkine(q);
Je = UR5e.jacob0(q);
% Je_dot = UR5e.jacob_dot(q,dq);
toc

% M = simplify(M);
% C = simplify(C);
% G = simplify(G);
Jb = simplify(Jb);
Je = simplify(Je);
g_st = simplify([g_st_.R, g_st_.t;
                 zeros(1,3),1]);
N = length(q);

Jb_dot = zeros(6,N);
Je_dot = zeros(6,N);

for k = 1 : N
    Jb_dot = Jb_dot + diff(Jb,q(k))*dq(k);
    Je_dot = Je_dot + diff(Je,q(k))*dq(k);
end

Jb_dot = simplify(Jb_dot);
Je_dot = simplify(Je_dot);
%%
% tic
% matlabFunction(M,'File','sub_direct/M_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'M_mat'});
% toc
% 
% tic
% matlabFunction(Jb,'File','sub_direct/Jb_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'J_mat'});
% toc
% 
% tic
% matlabFunction(Jb_dot,'File','sub_direct/Jb_dot_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6, dq1, dq2, dq3, dq4, dq5, dq6],'Outputs',{'dJ_mat'});
% toc

tic
matlabFunction(Je,'File','sub_direct/Je', ...
    'Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'Je_mat'});
toc

tic
matlabFunction(Je_dot,'File','sub_direct/Je_dot_fun', ...
    'Vars',[q1, q2, q3, q4, q5, q6, dq1, dq2, dq3, dq4, dq5, dq6],'Outputs',{'dJe_mat'});
toc



% tic
% matlabFunction(C,'File','sub_direct/C_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6, dq1, dq2, dq3, dq4, dq5, dq6],'Outputs',{'C_mat'});
% toc
% 
% tic
% matlabFunction(G,'File','sub_direct/G_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'G_torq'});
% toc
% 
% tic
% matlabFunction(g_st,'File','sub_direct/g_st_fun', ...
%     'Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'g_st_mat'});
% toc