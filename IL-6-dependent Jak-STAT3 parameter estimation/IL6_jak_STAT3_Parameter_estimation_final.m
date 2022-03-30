
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version of the Constrained Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
% Written by: ABHISEK BAKSHI
% Language: Matlab
% Date of creation: 06 May, 2019

%   - Number of States: 6
%   - Number of estimated prameters: 13
%   - Number of input to FIS: 2
%   - Number of output from FIS: 1
%
% works!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization

clear all;
close all;
clc;

% model parameters

x1	  =	0.02324  ;   % PG_3
x2	  =	0.07152	;	% glc_6P
x3	  =	0.06143	;	% glc arbitary
x4	  =	0.0922	;	% f_6P
x5	  =	0.06691	;	% f_16_BP 
x6	  =	0.01878	;	% f_26_BP arbitary



% Model used for estimation
model =@(t,y) IL6_jak_STAT3_extended_model(t,y);
jac=@(t,y) IL6_jak_STAT3_extended_model_jacobean(t,y);


% Set up time vector 
T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
time = [0:Ts:T]';


% Dimensions
N=19;


%%% Obtain measurement signal

% Output matrix

H = zeros(1,19);

H(1,3) = 1;
% H(2,1) = 1;

   
% Measurement noise covariance
R = eye(1)*0.001;
        
% Process noise
Q=0.01*eye(19);

% Initial conditions
xe0 = [0.07152	0.06143	0.0922	0.06691	0.01878	0.03506	0.0005922	0.0009618	0.0002324	0.0002419	0.00024725	0.0009615	0.0005368	0.0007354	0.0005286	0.0003492	0.0001416	0.000735	0.000631]';
Pe0=eye(19);


% State constraints Dx<=d
D = -eye(N);
d = -0.1*ones(N,1);

alpha = 0.56;%rand(1);

% Invoke CRFIEKF function
[estimated_parameter]= IL6_jak_STAT3_CRFIEKF_run(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha);
disp('Parameter estimation of IL6 dependent JAK/STAT pathway completed');
disp('Validation with previous model started... please wait...');
[x_v, reference] = IL6_jak_STAT3_simulation(estimated_parameter);
% disp('Simulation of JAK/STAT pathway using estimated parameter completed...');
disp('Generating plot...');
IL6_jak_STAT3_validation_plot(x_v, reference)



