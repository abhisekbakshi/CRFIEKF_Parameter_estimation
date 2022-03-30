function [xdot] = BPM_extended_model( t,x_v)

NOM = 11;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x1 = x_v(1);
x2 = x_v(2);
x3 = x_v(3);
x4 = x_v(4);


f13 = x_v(5);
f21 = x_v(6);
f32 = x_v(7);
f43 = x_v(8);
f44 = x_v(9);
f51 = x_v(10);
f64 = x_v(11);


g1 = 0.9;
g2 = 0.33;
g3 = 0.05;
g4 = 0.17;
g5 = 0.01;
g6 = 0.22;


x0 = 0.20;


%%%%% METABOLIC ODE PART

xdot(1) = (g1*x0)/(x3^f13) - g2*x1^f21 - g5*x1^f51;
xdot(2) = g2*x1^f21 - g3*x2^f32;
xdot(3) = g3*x2^f32 - g4*x3^f43*x4^f44;
xdot(4) = g5*x1^f51 - g6*x4^f64;
xdot(5) = 0;
xdot(6) = 0;
xdot(7) = 0;
xdot(8) = 0;
xdot(9) = 0;
xdot(10) = 0;
xdot(11) = 0;

end




