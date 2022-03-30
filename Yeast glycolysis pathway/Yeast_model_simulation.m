function [x_v, reference, reference_plosone_noiseless, reference_plosone_noisy] = Yeast_model_simulation(CRFIEKF_estimated_parameter)

warning('off','all');

global J0_CRFIEKF k1_CRFIEKF k2_CRFIEKF	k3_CRFIEKF k4_CRFIEKF k5_CRFIEKF k6_CRFIEKF	k_CRFIEKF kappa_CRFIEKF	q_CRFIEKF K1_CRFIEKF xi_CRFIEKF N_CRFIEKF A_CRFIEKF;

J0_CRFIEKF	=	CRFIEKF_estimated_parameter(1);
k1_CRFIEKF	=	CRFIEKF_estimated_parameter(2);
k2_CRFIEKF	=	CRFIEKF_estimated_parameter(3);
k3_CRFIEKF	=	CRFIEKF_estimated_parameter(4);
k4_CRFIEKF	=	CRFIEKF_estimated_parameter(5);
k5_CRFIEKF	=	CRFIEKF_estimated_parameter(6);
k6_CRFIEKF	=	CRFIEKF_estimated_parameter(7);
k_CRFIEKF	=	CRFIEKF_estimated_parameter(8);
kappa_CRFIEKF	=	CRFIEKF_estimated_parameter(9);
q_CRFIEKF	=	CRFIEKF_estimated_parameter(10);
K1_CRFIEKF	=	CRFIEKF_estimated_parameter(11);
xi_CRFIEKF	=	CRFIEKF_estimated_parameter(12);
N_CRFIEKF	=	CRFIEKF_estimated_parameter(13);
A_CRFIEKF	=	CRFIEKF_estimated_parameter(14);


x_v = zeros(1,7);
x_v(1:7) =   [0.501 1.955 0.198 0.148 0.161 0.161 0.064];


% uu = zeros(1, 27); 

T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';
tic

for ii = 1:(T/Ts)
                               
    [tsim, xsim] = ode23tb('Yeast_model', Tspan, x_v(ii,:));
    x_v(ii+1, :) = real(xsim(end, :));

end

disp('Yeast model simulation with parameter values by proposed model completed...');


reference = zeros(1,7);
reference(1:7) = [0.501 1.955 0.198 0.148 0.161 0.161 0.064];

for ii = 1:(T/Ts)
                               
    [tsim, xsim] = ode23tb('Yeast_model_original', Tspan, reference(ii,:));
    reference(ii+1, :) = real(xsim(end, :));

end

disp('Yeast model simulation with parameter values by Ruoff et al. (2003) completed...');


reference_plosone_noiseless = zeros(1,7);
reference_plosone_noiseless(1:7) = [0.501 1.955 0.198 0.148 0.161 0.161 0.064];

for ii = 1:(T/Ts)
                               
    [tsim, xsim] = ode23tb('Yeast_model_plosone_noiseless', Tspan, reference_plosone_noiseless(ii,:));
    reference_plosone_noiseless(ii+1, :) = real(xsim(end, :));

end

disp('Yeast model simulation with parameter values by Yazdani et al. (2020) noiseless data completed...');




reference_plosone_noisy = zeros(1,7);
reference_plosone_noisy(1:7) = [0.501 1.955 0.198 0.148 0.161 0.161 0.064];

for ii = 1:(T/Ts)
                               
    [tsim, xsim] = ode23tb('Yeast_model_plosone_noisy', Tspan, reference_plosone_noisy(ii,:));
    reference_plosone_noisy(ii+1, :) = real(xsim(end, :));

end

disp('Yeast model simulation with parameter values by Yazdani et al. (2020) noisy data completed...');


end


