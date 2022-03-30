function [xdot] = BPM_model_original(t,x_v)

NOM = 4;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x0 = 0.6;

x1 = x_v(1);
x2 = x_v(2);
x3 = x_v(3);
x4 = x_v(4);


f13 = 0.8;
f21 = 0.5;
f32 = 0.75;
f43 = 0.5;
f44 = 0.2;
f51 = 0.4999;
f64 = 0.7999;


g1 = 19.9999;
g2 = 7.9998;
g3 = 2.9998;
g4 = 4.9998;
g5 = 2.0002;
g6 = 5.9997;




%%%%% METABOLIC ODE PART

xdot(1) = (g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51;
xdot(2) = g2*x1^f21 - g3*x2^f32;
xdot(3) = g3*x2^f32 - g4*x3^f43*x4^f44;
xdot(4) = g5*x1^f51 - g6*x4^f64;




end