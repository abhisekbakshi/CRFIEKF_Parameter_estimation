function [xdot] = Yeast_model_plosone_noiseless( t,x_v)

NOM = 7;              % the number of molecules
xdot = zeros(NOM, 1);


S1	=	x_v(1);     				
S2	=	x_v(2);    	
S3	=	x_v(3);     				
S4	=	x_v(4);     				
S5	=	x_v(5);    			
S6	=	x_v(6);     				
S7	=	x_v(7);    

J0	=	0.02409759	;
k1	=	0.9980002	;
k2	=	0.059194081	;
k3	=	0.158084192	;
k4	=	1	;
k5	=	0.01189881	;
k6	=	0.119088091	;
k	=	0.0169983	;
kappa	=	0.129087091	;
q	=	0.03909609	;
K1	=	0.00429957	;
xi	=	9.39906E-05	;
N	=	0.009089091	;
A	=	0.03919608	;



%%%%% METABOLIC ODE PART


xdot(1) = J0 - (k1*S1*S6)/(1 + S6/K1)^q;
xdot(2) = 2*(k1*S1*S6)/(1 + S6/k1)^q - k2*S2*(N-S5) - k6*S2*S5;
xdot(3) = k2*S2*(N-S5) - k3*S3*(A - S6);
xdot(4) = k3*S3*(A - S6) - k4*S4*S5 -kappa*(S4 - S7);
xdot(5) = k2*S2*(N-S5) - k4*S4*S5 - k6*S2*S5;
xdot(6) = -2*(k1*S1*S6)/(1 + S6/k1)^q + 2*k3*S3*(A - S6) - k5*S6;
xdot(7) = xi*k*(S4 - S7) - k*S7;

end






