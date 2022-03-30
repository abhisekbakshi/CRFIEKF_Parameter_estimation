function jac = Yeast_model_jacobean(t,x_v)


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


% Preallocate Jacobian
jac=zeros(21,21);
% Set nonzero elements



jac(1, 1) = -(S6*k1)/(S6/K1 + 1)^q;

jac(1, 6) = -(S1*k1*(K1 + S6 - S6*q))/((S6/K1 + 1)^q*(K1 + S6));

jac(1, 8) = 1;

jac(1, 9) = -(S1*S6)/(S6/K1 + 1)^q;

jac(1, 17) = (S1*S6*k1*log(S6/K1 + 1))/(S6/K1 + 1)^q;

jac(1, 18) = -(S1*S6^2*k1*q)/(K1^2*(S6/K1 + 1)^(q + 1));

jac(2, 1) = (2*S6*k1)/(S6/k1 + 1)^q;

jac(2, 2) = - S5*k6 - k2*(N - S5);

jac(2, 5) = S2*(k2 - k6);

jac(2, 6) = (2*S1*k1*(S6 + k1 - S6*q))/((S6 + k1)*(S6/k1 + 1)^q);

jac(2, 9) = (2*S1*S6*(S6 + k1 + S6*q))/((S6 + k1)*(S6/k1 + 1)^q);

jac(2, 10) = -S2*(N - S5);

jac(2, 14) = -S2*S5;

jac(2, 17) = -(2*S1*S6*k1*log(S6/k1 + 1))/(S6/k1 + 1)^q;

jac(2, 20) = -S2*k2;

jac(3, 2) = k2*(N - S5);

jac(3, 3) = -k3*(A - S6);

jac(3, 5) = -S2*k2;

jac(3, 6) = S3*k3;

jac(3, 10) = S2*(N - S5);

jac(3, 11) = -S3*(A - S6);

jac(3, 20) = S2*k2;

jac(3, 21) = -S3*k3;

jac(4, 3) = k3*(A - S6);

jac(4, 4) = - kappa - S5*k4;

jac(4, 5) = -S4*k4;

jac(4, 6) = -S3*k3;

jac(4, 7) = kappa;

jac(4, 11) = S3*(A - S6);

jac(4, 12) = -S4*S5;

jac(4, 16) = S7 - S4;

jac(4, 21) = S3*k3;

jac(5, 2) = k2*(N - S5) - S5*k6;

jac(5, 4) = -S5*k4;

jac(5, 5) = - S2*k2 - S2*k6 - S4*k4;

jac(5, 10) = S2*(N - S5);

jac(5, 12) = -S4*S5;

jac(5, 14) = -S2*S5;

jac(5, 20) = S2*k2;

jac(6, 1) = -(2*S6*k1)/(S6/k1 + 1)^q;

jac(6, 3) = 2*k3*(A - S6);

jac(6, 6) = (2*S1*S6*q)/(S6/k1 + 1)^(q + 1) - 2*S3*k3 - (2*S1*k1)/(S6/k1 + 1)^q - k5;

jac(6, 9) = -(2*S1*S6*(S6 + k1 + S6*q))/((S6 + k1)*(S6/k1 + 1)^q);

jac(6, 11) = 2*S3*(A - S6);

jac(6, 13) = -S6;

jac(6, 17) = (2*S1*S6*k1*log(S6/k1 + 1))/(S6/k1 + 1)^q;

jac(6, 21) = 2*S3*k3;

jac(7, 4) = k*xi;

jac(7, 7) = -k*(xi + 1);

jac(7, 15) = xi*(S4 - S7) - S7;

jac(7, 19) = k*(S4 - S7);


end



