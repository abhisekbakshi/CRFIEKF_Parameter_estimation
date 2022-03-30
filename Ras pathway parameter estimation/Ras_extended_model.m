function [xdot] = Ras_extended_model( t,x_v)

NOM = 22;              % the number of molecules
xdot = zeros(NOM, 1);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

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
k1  = x_v(12);
k2  = x_v(13);
k3  = x_v(14);
k4  = x_v(15);
k5  = x_v(16);
k6  = x_v(17);
k7  = x_v(18);
k8  = x_v(19);
k9  = x_v(20);
k10 = x_v(21);
k11 = x_v(22);



%%%%% SIGNALING ODE PART

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
xdot(12) = 0;
xdot(13) = 0;
xdot(14) = 0;
xdot(15) = 0;
xdot(16) = 0;
xdot(17) = 0;
xdot(18) = 0;
xdot(19) = 0;
xdot(20) = 0;
xdot(21) = 0;
xdot(22) = 0;

 
