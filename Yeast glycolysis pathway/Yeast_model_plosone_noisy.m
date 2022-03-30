function [xdot] = Yeast_model_plosone_noisy( t,x_v)

NOM = 7;              % the number of molecules
xdot = zeros(NOM, 1);


S1	=	x_v(1);     				
S2	=	x_v(2);    	
S3	=	x_v(3);     				
S4	=	x_v(4);     				
S5	=	x_v(5);    			
S6	=	x_v(6);     				
S7	=	x_v(7);    

J0	=	0.024611928	;
k1	=	0.878472222	;
k2	=	0.04564951	;
k3	=	0.142156863	;
k4	=	0.990808824	;
k5	=	0.011846405	;
k6	=	0.128880719	;
k	=	0.015012255	;
kappa	=	0.136029412	;
q	=	0.040747549	;
K1	=	0.004799837	;
xi	=	2.34886E-05	;
N	=	0.012357026	;
A	=	0.042585784	;




%%%%% METABOLIC ODE PART


xdot(1) = J0 - (k1*S1*S6)/(1 + S6/K1)^q;
xdot(2) = 2*(k1*S1*S6)/(1 + S6/k1)^q - k2*S2*(N-S5) - k6*S2*S5;
xdot(3) = k2*S2*(N-S5) - k3*S3*(A - S6);
xdot(4) = k3*S3*(A - S6) - k4*S4*S5 -kappa*(S4 - S7);
xdot(5) = k2*S2*(N-S5) - k4*S4*S5 - k6*S2*S5;
xdot(6) = -2*(k1*S1*S6)/(1 + S6/k1)^q + 2*k3*S3*(A - S6) - k5*S6;
xdot(7) = xi*k*(S4 - S7) - k*S7;

end






