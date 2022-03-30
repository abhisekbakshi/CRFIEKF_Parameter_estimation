function [xdot] = JAK_STAT_extended_model( t,x_v)

NOM = 8;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x1 = x_v(1);
x2 = x_v(2);
x3 = x_v(3);
x4 = x_v(4);

k1 = x_v(5);
k2 = x_v(6);
k3 = x_v(7);
k4 = x_v(8);



u = (2*k4*x3)/(k1*x1);

%%%%% METABOLIC ODE PART

xdot(1) = 2*k4*x3 - k1*u*x1;
xdot(2) = k1*u*x1 - 2*k2*x2^2;
xdot(3) = k2*x2^2 - k3*x3;
xdot(4) = k3*x3 - k4*x3;
xdot(5) = 0;
xdot(6) = 0;
xdot(7) = 0;
xdot(8) = 0;