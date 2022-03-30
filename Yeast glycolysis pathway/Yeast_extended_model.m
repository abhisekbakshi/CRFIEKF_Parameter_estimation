function [xdot] = Yeast_extended_model( t,x_v)

NOM = 21;              % the number of molecules
xdot = zeros(NOM, 1);

S1	=	x_v(1);     				
S2	=	x_v(2);    	
S3	=	x_v(3);     				
S4	=	x_v(4);     				
S5	=	x_v(5);    			
S6	=	x_v(6);     				
S7	=	x_v(7);    

J0	=	x_v(8);     					
k1	=	x_v(9);	    					
k2	=	x_v(10);						
k3	=	x_v(11);					
k4	=	x_v(12);					
k5	=	x_v(13);			
k6	=	x_v(14);					
k	=	x_v(15);			
kappa	=	x_v(16);						
q	=	x_v(17);		
K1	=	x_v(18);	
xi	=	x_v(19);					
N	=	x_v(20);				
A	=	x_v(21);					




%%%%% METABOLIC ODE PART


xdot(1) = J0 - (k1*S1*S6)/(1 + S6/K1)^q;
xdot(2) = 2*(k1*S1*S6)/(1 + S6/k1)^q - k2*S2*(N-S5) - k6*S2*S5;
xdot(3) = k2*S2*(N-S5) - k3*S3*(A - S6);
xdot(4) = k3*S3*(A - S6) - k4*S4*S5 -kappa*(S4 - S7);
xdot(5) = k2*S2*(N-S5) - k4*S4*S5 - k6*S2*S5;
xdot(6) = -2*(k1*S1*S6)/(1 + S6/k1)^q + 2*k3*S3*(A - S6) - k5*S6;
xdot(7) = xi*k*(S4 - S7) - k*S7;
xdot(8) = 0;
xdot(9) = 0;
xdot(10) = 0;
xdot(11) = 0;
xdot(12) = 0;
xdot(13) = 0;
xdot(14) = 0;
xdot(15) = 0;
xdot(16) = 0;
xdot(17) = 0;
xdot(18) = 0;
xdot(19) = 0;
xdot(21) = 0;


end


