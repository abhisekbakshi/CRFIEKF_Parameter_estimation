%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version of the Constrained Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
% Written by: ABHISEK BAKSHI
% Language: Matlab
% Date of creation: 07 Dec, 2021

%   - Number of States: 11
%   - Number of estimated prameters: 11
%   - Number of input to FIS: 2
%   - Number of output from FIS: 2
%
% works!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Initialization

clear all;
close all;
clc;

warning('off','all');

% model parameters
x1 = 66;
x2 = 0.054;
x3 = 0.019;
x4 = 59;
x5 = 0.09;
x6 = 0.012;
x7 = 65;
x8 = 26;
x9 = 175;
x10 = 161;
x11 = 2.18;
k1 = 0.546;
k2 = 0.014;
k3 = 0.619;
k4 = 0.046;
k5 = -1.29;
k6 = 0.84;
k7 = -0.05;
k8 = 0.43;
k9 = 0.98;
k10 = -0.006;
k11 = 0.88;



p=[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11];

%%% control input EpoR
% u = [0,9,14,44,59,50,46,32,35,19,19,17,3,2,1,1];


% Model used for estimation
model =@(t,y) Ras_extended_model(t,y);
jac=@(t,y) Ras_extended_model_jacobean(t,y);



% % % % % Find the initial conditions (equilibrium point)
x0 = [66	0.054	0.019	59	0.09	0.012	65	26	175	 161	 2.18];

% Dimensions
S=length(x0);
P=length(p);
N=S+P;


% Output matrix
H=[0 0 0 0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 0 0 1 0 0];

H=[H zeros(size(H,1),P)];

time = [0 1 2 3 4 5 6 7 8 9 10];
        
% Measurement noise covariance
R = eye(2)*0.01;

% Process noise
Q=0.01*eye(22);

% Initial conditions
xe0 = [66	0.054	0.019	59	0.09	0.012	65	26	175	 161	 2.18  0.53	0.0072	0.625	0.00105	0.0100	0.8	0.35	0.071	0.92	0.00122	0.87];
Pe0=eye(22);

alpha = 0.1;
reg_mat = Pe0;

[estimated_parameter]=CRFIEKF_Ras (time, model, jac, H , Q, R, xe0, Pe0, alpha);
disp('Parameter estimation of Ras pathway completed ')
disp('Validation with previous result started... please wait...');
[x_v] = Ras_model_simulation(estimated_parameter);
disp('Validation with WET lab data completed');
disp('Generating plot...');
Ras_model_validation_plot(x_v);

