function [xdot] = JAK_STAT_model( t,x_v)

NOM = 5;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%
%%% m39 is not included intentionally

x1 = x_v(1);
x2 = x_v(2);
x3 = x_v(3);
x4 = x_v(4);
u =  x_v(5);

global k1 k2 k3 k4;

k1 = 0.0201;


%%%%% METABOLIC ODE PART

xdot(1) = 2*k4*x4 - k1*u*x1;
xdot(2) = k1*u*x1 - 2*k2*x2^2;
xdot(3) = k2*x2^2 - k3*x3;
xdot(4) = k3*x3 - k4*x4;
xdot(5) = 0;


end








