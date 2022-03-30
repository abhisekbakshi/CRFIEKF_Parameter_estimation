
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version of the Constrained Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
% Written by: ABHISEK BAKSHI
% Language: Matlab
% Date of creation: 06 May, 2019

%   - Number of States: 4
%   - Number of estimated prameters: 7
%   - Number of input to FIS: 3
%   - Number of output from FIS: 2
%
% works!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

warning('off','all');

% model parameters
x0 = 0.20;
x1 = 0.47;
x2 = 0.9;
x3 = 0.4;
x4 = 0.13;

g1 = 0.9;
g2 = 0.33;
g3 = 0.05;
g4 = 0.17;
g5 = 0.01;
g6 = 0.22;



f13 = 0.2;
f21 = 0.6;
f32 = 0.1;
f43 = 0.1;
f44 = 0.1;
f51 = 0.1;
f64 = 0.1;



% Model used for estimation
model =@(t,x_v) BPM_extended_model(t,x_v);
jac=@(t,x_v) BPM_extended_model_jacobean(t,x_v);

N = 11;

% Output matrix
H = zeros(2,11);

H(1,1) = 1;
H(2,3) = 1;


% Set up time vector 
T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
time = [0:Ts:T]';


% Measurement noise covariance
R = eye(1)*0.001;
        
% Process noise
Q=0.1*eye(11);

% Initial conditions
xe0 = [x1, x2, x3, x4, f13, f21, f32, f43, f44, f51, f64];
Pe0=eye(11);

% State constraints Dx<=d
D = -eye(N);
d = -0.01*ones(N,1);

alpha = 0.56;

[estimated_parameter]=CRFIEKF_BPM(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha, T, Ts);
disp('Parameter estimation completed');
disp('Validation with previous model started... please wait...');
[x_v, reference] = BPM_model_simulation(estimated_parameter);
disp('Validation with original model completed');
disp('Generating plot...');
BPM_model_validation_plot(x_v, reference);


