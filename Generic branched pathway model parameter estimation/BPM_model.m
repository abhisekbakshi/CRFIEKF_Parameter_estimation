function [xdot] = BPM_model(t,x_v)

NOM = 4;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x0 = 0.6;

x1 = x_v(1);
x2 = x_v(2);
x3 = x_v(3);
x4 = x_v(4);


global f13 f21 f32 f43 f44 f51 f64;

g1 = 20;
g2 = 8;
g3 = 3;
g4 = 5;
g5 = 2;
g6 = 6;


f64 = 0.8;


%%%%% METABOLIC ODE PART

xdot(1) = (g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51;
xdot(2) = g2*x1^f21 - g3*x2^f32;
xdot(3) = g3*x2^f32 - g4*x3^f43*x4^f44;
xdot(4) = g5*x1^f51 - g6*x4^f64;




end