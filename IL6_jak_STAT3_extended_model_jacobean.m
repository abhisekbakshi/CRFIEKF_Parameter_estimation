function jac = IL6_jak_STAT3_extended_model_jacobean_final(t,x_v)


x1	=	x_v(1);     %	IL-6IL-6R
x2	=	x_v(2);     %	 gp130
x3	=	x_v(3);     %	(p)Rcomplex
x4	=	x_v(4);     %	(p)STAT3
x5	=	x_v(5);     %	SOCS3 mRNA
x6	=	x_v(6);     %	SOCS3 protein

p1	=	x_v(7);     
p2	=	x_v(8);     
p3	=	x_v(9);     
p4	=	x_v(10);     
p5	=	x_v(11);     
p6	=	x_v(12);     
p7	=	x_v(13);     
p8	=	x_v(14);     
p9	=	x_v(15);    
p10	=	x_v(16);     
p11	=	x_v(17);     
p12	=	x_v(18);     
p13	=	x_v(19);     

u = 0.17;


x8 = (16.8 - x2 - 2*x3)/2;      %  Rcomplex
x7 = 2.1 - x1 - 2*x3 - 2*x8;    %  IL-6R
x9 = 83 - x4;                   %  STAT3

%%%%% ENZYME/GENE INITIAL CONCENTRATION g in nM  %%%%%

% Preallocate Jacobian
jac=zeros(19,19);
% Set nonzero elements



jac(1, 1) = - p2 - 4*p3*x1*x2^2;

jac(1, 2) = -4*p3*x1^2*x2;

jac(1, 7) = u*x7;

jac(1, 8) = -x1;

jac(1, 9) = -2*x1^2*x2^2;

jac(1, 10) = 2*x8;

jac(2, 1) = -4*p3*x1*x2^2;

jac(2, 2) = -4*p3*x1^2*x2;

jac(2, 9) = -2*x1^2*x2^2;

jac(2, 10) = 2*x8;

jac(3, 3) = -p6;

jac(3, 6) = -(p5*p13*x8)/(p13*x6 + 1)^2;

jac(3, 11) = x8/(p13*x6 + 1);

jac(3, 12) = -x3;

jac(3, 19) = -(p5*x6*x8)/(p13*x6 + 1)^2;

jac(4, 3) = p7*x9;

jac(4, 4) = -p8;

jac(4, 13) = x3*x9;

jac(4, 14) = -x4;

jac(5, 4) = p9;

jac(5, 5) = -p10;

jac(5, 15) = x4;

jac(5, 16) = -x5;

jac(6, 5) = p11;

jac(6, 6) = -p12;

jac(6, 17) = x5;

jac(6, 18) = -x6;



