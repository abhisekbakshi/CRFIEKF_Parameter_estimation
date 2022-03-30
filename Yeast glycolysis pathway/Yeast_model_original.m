function [xdot] = Yeast_model_original( t,x_v)

NOM = 7;              % the number of molecules
xdot = zeros(NOM, 1);


S1	=	x_v(1);     				
S2	=	x_v(2);    	
S3	=	x_v(3);     				
S4	=	x_v(4);     				
S5	=	x_v(5);    			
S6	=	x_v(6);     				
S7	=	x_v(7); 


J0	=	0.024512256	;
k1	=	1	;
k2	=	0.059529765	;
k3	=	0.15957979	;
k4	=	1	;
k5	=	0.012306153	;
k6	=	0.11955978	;
k	=	0.017508754	;
kappa	=	0.129564782	;
q	=	0.03951976	;
K1	=	0.004702351	;
xi	=	0.00050025	;
N	=	0.009504752	;
A	=	0.03951976	;


%%%%% METABOLIC ODE PART


xdot(1) = J0 - (k1*S1*S6)/(1 + S6/K1)^q;
xdot(2) = 2*(k1*S1*S6)/(1 + S6/k1)^q - k2*S2*(N-S5) - k6*S2*S5;
xdot(3) = k2*S2*(N-S5) - k3*S3*(A - S6);
xdot(4) = k3*S3*(A - S6) - k4*S4*S5 -kappa*(S4 - S7);
xdot(5) = k2*S2*(N-S5) - k4*S4*S5 - k6*S2*S5;
xdot(6) = -2*(k1*S1*S6)/(1 + S6/k1)^q + 2*k3*S3*(A - S6) - k5*S6;
xdot(7) = xi*k*(S4 - S7) - k*S7;

end






