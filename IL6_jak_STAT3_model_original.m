function [xdot] = IL6_jak_STAT3_model_original( t,x_v)

NOM = 6;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x1	=	x_v(1);     %	IL-6IL-6R
x2	=	x_v(2);     %	 gp130
x3	=	x_v(3);     %	(p)Rcomplex
x4	=	x_v(4);     %	(p)STAT3
x5	=	x_v(5);     %	SOCS3 mRNA
x6	=	x_v(6);     %	SOCS3 protein

global x8;
global x7;
global x9;

x8 = (16.8 - x2 - 2*x3)/2;      %  Rcomplex
x7 = 2.1 - x1 - 2*x3 - 2*x8;    %  IL-6R
x9 = 83 - x4;                   %  STAT3

u = 0.17;

%%%%% Kinetic parameter CONSTANT km in nM %%%%%


p1 =    0.1218;
p2 =    0.0388;
p3 =    3.5920;
p4 =    0.0484;
p5 =    0.0803;
p6 =    0.0864;
p7 =    0.156;
p8 =    0.010;
p9 =    0.026;
p10 =   0.021;
p11 =   0.029;
p12 =   0.008;
p13 =   0.029;



%%%%% METABOLIC ODE PART


xdot(1) = p1*x7*u - p2*x1 - 2*p3*x2^2*x1^2 + 2*p4*x8;
xdot(2) = 2*p4*x8 - 2*p3*x2^2*x1^2;
xdot(3) = (p5*x8)/(1 + p13*x6) - p6*x3;
xdot(4) = p7*x3*x9 - p8*x4;
xdot(5) = p9*x4 - p10*x5;
xdot(6) = p11*x5 - p12*x6;

end








