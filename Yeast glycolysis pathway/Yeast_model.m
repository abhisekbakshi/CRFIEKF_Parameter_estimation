function [xdot] = Yeast_model( t,x_v)

NOM = 7;              % the number of molecules
xdot = zeros(NOM, 1);


S1	=	x_v(1);     				
S2	=	x_v(2);    	
S3	=	x_v(3);     				
S4	=	x_v(4);     				
S5	=	x_v(5);    			
S6	=	x_v(6);     				
S7	=	x_v(7);    


global J0_CRFIEKF k1_CRFIEKF k2_CRFIEKF	k3_CRFIEKF k4_CRFIEKF k5_CRFIEKF k6_CRFIEKF	k_CRFIEKF kappa_CRFIEKF	q_CRFIEKF K1_CRFIEKF xi_CRFIEKF N_CRFIEKF A_CRFIEKF;


q_CRFIEKF	=	0.187972548	;
N_CRFIEKF	=	0.0427431414	;


%%%%% METABOLIC ODE PART


xdot(1) = J0_CRFIEKF - (k1_CRFIEKF*S1*S6)/(1 + S6/K1_CRFIEKF)^q_CRFIEKF;
xdot(2) = 2*(k1_CRFIEKF*S1*S6)/(1 + S6/k1_CRFIEKF)^q_CRFIEKF - k2_CRFIEKF*S2*(N_CRFIEKF-S5) - k6_CRFIEKF*S2*S5;
xdot(3) = k2_CRFIEKF*S2*(N_CRFIEKF-S5) - k3_CRFIEKF*S3*(A_CRFIEKF - S6);
xdot(4) = k3_CRFIEKF*S3*(A_CRFIEKF - S6) - k4_CRFIEKF*S4*S5 -kappa_CRFIEKF*(S4 - S7);
xdot(5) = k2_CRFIEKF*S2*(N_CRFIEKF-S5) - k4_CRFIEKF*S4*S5 - k6_CRFIEKF*S2*S5;
xdot(6) = -2*(k1_CRFIEKF*S1*S6)/(1 + S6/k1_CRFIEKF)^q_CRFIEKF + 2*k3_CRFIEKF*S3*(A_CRFIEKF - S6) - k5_CRFIEKF*S6;
xdot(7) = xi_CRFIEKF*k_CRFIEKF*(S4 - S7) - k_CRFIEKF*S7;

end






