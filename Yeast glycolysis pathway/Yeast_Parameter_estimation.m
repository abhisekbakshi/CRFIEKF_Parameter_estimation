
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version of the Constrained Regularized Fuzzy Inferred Extended Kalman Filter for parameter estimation of glycolysis pathway
% Written by: ABHISEK BAKSHI
% Language: Matlab
% Date of creation: 06 May, 2019

%   - Number of States: 7
%   - Number of estimated prameters: 14
%   - Number of input to FIS: 4
%   - Number of output from FIS: 2
%
% works!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization

clear all;
close all;
clc;

% model parameters

S1	=	0.3501;     				
S2	=	0.955;    	
S3	=	0.0198;     				
S4	=	0.0148;     				
S5	=	0.0161;    			
S6	=	0.0161;     				
S7	=	0.0064; 


J0	=	0.0124512256	;
k1	=	0.9	;
k2	=	0.0459529765	;
k3	=	0.105957979	;
k4	=	0.9	;
k5	=	0.00912306153	;
k6	=	0.0911955978	;
k	=	0.00917508754	;
kappa	=	0.09129564782	;
q	=	0.023951976	;
K1	=	0.0034702351	;
xi	=	0.000450025	;
N	=	0.0089504752	;
A	=	0.023951976	;


%%% Number of parameters = 73 
p = [S1, S2, S3, S4, S5, S6, S7, J0, k1, k2, k3, k4, k5, k6, k, kappa, q, K1, xi, N, A];

% Model used for estimation
model =@(t,y) Yeast_extended_model(t,y);
jac=@(t,y) Yeast_model_jacobean(t,y);


% Set up time vector 
T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
time = [0:Ts:T]';

% Dimensions
Nd = 21;


%%% Obtain measurement signal

% Output matrix

H = zeros(2,21);

H(1,2) = 1;
H(2,4) = 1;


% Measurement noise covariance
R = eye(2)*0.001;
        
% Process noise
Q=0.01*eye(21);

% Initial conditions
xe0 = [S1, S2, S3, S4, S5, S6, S7, J0, k1, k2, k3, k4, k5, k6, k, kappa, q, K1, xi, N, A]';
Pe0=eye(21);


% State constraints Dx<=d
D = -eye(Nd);
d = -0.1*ones(Nd,1);

alpha = 0.56;%rand(1);

% Invoke CRFIEKF function
[Estimated_parameter]= Yeast_glycolysis_CRFIEKF(time, model, jac, H , Q, R, xe0, Pe0, D, d, alpha);
disp('Yeast model glycolysis Parameter estimation completed');
disp('Validation with Ruoff et al. (2003) and Yazdani et al. (2020) started... please wait...');
[x_v, reference, reference_plosone_noiseless, reference_plosone_noisy] = Yeast_model_simulation(Estimated_parameter);
disp('Validation with Ruoff et al. (2003) and (Yazdani et al. (2020) completed');
disp('Generating plot...');
Yeast_glyolysis_validation_plot(x_v, reference, reference_plosone_noiseless, reference_plosone_noisy);


