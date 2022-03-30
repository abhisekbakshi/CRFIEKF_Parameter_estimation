function [xdot] = Ras_model_for_validation( t,x_v)

NOM = 11;         % the number of molecules
xdot = zeros(NOM, 1);
global k1 k2 k3 k4 k5 k6 k7 k8 k9 k10  k11;

x1  = x_v(1);
x2  = x_v(2);
x3  = x_v(3);
x4  = x_v(4);
x5  = x_v(5);
x6  = x_v(6);
x7  = x_v(7);
x8  = x_v(8);
x9  = x_v(9);
x10 = x_v(10);
x11 = x_v(11);


%%%%% METABOLIC ODE PART

xdot(1) = k2*x3 + k5*x4 - k1*x1*x2;
xdot(2) = k2*x3 + k11*x11 - k1*x1*x2;
xdot(3) = k4*x4 - k2*x3 + k1*x1*x2 - k3*x3*x9;
xdot(4) = k3*x3*x9 - k5*x4 - k4*x4;
xdot(5) = k5*x4 + k7*x8 - k6*x5*x7;
xdot(6) = k5*x4 + k10*x11 - k9*x6*x10;
xdot(7) = k7*x8 + k8*x8 - k6*x5*x7;
xdot(8) = k6*x5*x7 - k8*x8 - k7*x8;
xdot(9) = k4*x4 + k8*x8 - k3*x3*x9;
xdot(10) = k10*x11 + k11*x11 - k9*x6*x10;
xdot(11) = k9*x6*x10 - k11*x11 - k10*x11;
