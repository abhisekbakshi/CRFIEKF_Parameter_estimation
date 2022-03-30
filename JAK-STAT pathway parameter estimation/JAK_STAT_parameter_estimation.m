%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version of the Constrained Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
% Written by: ABHISEK BAKSHI
% Language: Matlab
% Date of creation: 06 May, 2019

%   - Number of States: 4
%   - Number of estimated prameters: 4
%   - Number of input to FIS: 1
%   - Number of output from FIS: 1
%
% works!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Initialization

clear all;
close all;
clc;

warning('off','all');

% model parameters
x1 = 0.1;
x2 = 0;
x3 = 0;
x4 = 0;


k1 = 0.017;
k2 = 1.1768;
k3 = 0.1184;
k4 = 0.1;
p=[k1 k2 k3 k4];

%%% control input EpoR
u = [0,9,14,44,59,50,46,32,35,19,19,17,3,2,1,1];

% Model used for estimation
model =@(t,y) JAK_STAT_extended_model(t,y);
jac=@(t,y) JAK_STAT_model_jacobian(t,y);

% % Initial conditions (equilibrium point)
x0 = [0.1	0	0	0];

% Dimensions
S=length(x0);
P=length(p);
N=S+P;

% Output matrix
H=[1 0 0 0];
H=[H zeros(size(H,1),P)];


% Set up time vector 
T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
time = [0:Ts:T]';

% time  = [0,2,4,6,8,10,12,14,16,18,20,25,30,40,50,60];
        
% Measurement noise covariance
R = eye(1)*0.001;
        
% Process noise
Q=0.1*eye(8);

% Initial conditions
xe0 = [0.1	0	0	0	0.017	2.2	0.1184	0.1];
Pe0=eye(8);

% State constraints Dx<=d
D = -eye(N);
d = -0.01*ones(N,1);

alpha = 0.15;

[estimated_parameter]= CRFIEKF_JAK_STAT (time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha);
disp('Parameter estimation of JAK/STAT pathway completed...');
pause(5);
[x_v, ph_stat5, total_stat5] = JAK_STAT_simulation(estimated_parameter);
total_stat5 = total_stat5';
ph_stat5 = ph_stat5';
disp('Simulation of JAK/STAT pathway using estimated parameter completed...');
disp('Generating plot...');
pause(3);
JA_STAT_model_validation_plot(x_v, total_stat5, ph_stat5);






