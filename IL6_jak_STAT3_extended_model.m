function [xdot] = IL6_jak_STAT3_extended_model( t,x_v)

NOM = 19;              % the number of molecules
xdot = zeros(NOM, 1);


x1	=	x_v(1);     %	IL-6IL-6R
x2	=	x_v(2);     %	 gp130
x3	=	x_v(3);     %	(p)Rcomplex
x4	=	x_v(4);     %	(p)STAT3
x5	=	x_v(5);     %	SOCS3 mRNA
x6	=	x_v(6);     %	SOCS3 protein

p1	=	x_v(7);     
p2	=	x_v(8);     
p3	=	3.8920; 
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

x8 = (16.8 - x2 - 2*x3)/2;      %  Rcomplex
x7 = 2.1 - x1 - 2*x3 - 2*x8;    %  IL-6R
x9 = 83 - x4;                   %  STAT3
u=0.17;

%%%%% METABOLIC ODE PART


xdot(1) = p1*x7*u - p2*x1 - 2*p3*x2^2*x1^2 + 2*p4*x8;
xdot(2) = 2*p4*x8 - 2*p3*x2^2*x1^2;
xdot(3) = (p5*x8)/(1 + p13*x6) - p6*x3;
xdot(4) = p7*x3*x9 - p8*x4;
xdot(5) = p9*x4 - p10*x5;
xdot(6) = p11*x5 - p12*x6;
xdot(7) = 0;
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

end

