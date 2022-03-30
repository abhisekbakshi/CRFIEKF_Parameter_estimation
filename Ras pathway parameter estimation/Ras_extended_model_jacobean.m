function jac = Ras_extended_model_jacobean(t,y)


x1  = y(1);
x2  = y(2);
x3  = y(3);
x4  = y(4);
x5  = y(5);
x6  = y(6);
x7  = y(7);
x8  = y(8);
x9  = y(9);
x10 = y(10);
x11 = y(11);
k1  = y(12);
k2  = y(13);
k3  = y(14);
k4  = y(15);
k5  = y(16);
k6  = y(17);
k7  = y(18);
k8  = y(19);
k9  = y(20);
k10 = y(21);
k11 = y(22);


% Preallocate Jacobian
jac=zeros(22,22);
% Set nonzero elements

jac(1, 1) = -k1*x2;

jac(1, 2) = -k1*x1;

jac(1, 3) = k2;

jac(1, 4) = k5;

jac(1, 12) = -x1*x2;

jac(1, 13) = x3;

jac(1, 16) = x4;

jac(2, 1) = -k1*x2;

jac(2, 2) = -k1*x1;

jac(2, 3) = k2;

jac(2, 11) = k11;

jac(2, 12) = -x1*x2;

jac(2, 13) = x3;

jac(2, 22) = x11;

jac(3, 1) = k1*x2;

jac(3, 2) = k1*x1;

jac(3, 3) = - k2 - k3*x9;

jac(3, 4) = k4;

jac(3, 9) = -k3*x3;

jac(3, 12) = x1*x2;

jac(3, 13) = -x3;

jac(3, 14) = -x3*x9;

jac(3, 15) = x4;

jac(4, 3) = k3*x9;

jac(4, 4) = - k4 - k5;

jac(4, 9) = k3*x3;

jac(4, 14) = x3*x9;

jac(4, 15) = -x4;

jac(4, 16) = -x4;

jac(5, 4) = k5;

jac(5, 5) = -k6*x7;

jac(5, 7) = -k6*x5;

jac(5, 8) = k7;

jac(5, 16) = x4;

jac(5, 17) = -x5*x7;

jac(5, 18) = x8;

jac(6, 4) = k5;

jac(6, 6) = -k9*x10;

jac(6, 10) = -k9*x6;

jac(6, 11) = k10;

jac(6, 16) = x4;

jac(6, 20) = -x6*x10;

jac(6, 21) = x11;

jac(7, 5) = -k6*x7;

jac(7, 7) = -k6*x5;

jac(7, 8) = k7 + k8;

jac(7, 17) = -x5*x7;

jac(7, 18) = x8;

jac(7, 19) = x8;

jac(8, 5) = k6*x7;

jac(8, 7) = k6*x5;

jac(8, 8) = - k7 - k8;

jac(8, 17) = x5*x7;

jac(8, 18) = -x8;

jac(8, 19) = -x8;

jac(9, 3) = -k3*x9;

jac(9, 4) = k4;

jac(9, 8) = k8;

jac(9, 9) = -k3*x3;

jac(9, 14) = -x3*x9;

jac(9, 15) = x4;

jac(9, 19) = x8;

jac(10, 6) = -k9*x10;

jac(10, 10) = -k9*x6;

jac(10, 11) = k10 + k11;

jac(10, 20) = -x6*x10;

jac(10, 21) = x11;

jac(10, 22) = x11;

jac(11, 6) = k9*x10;

jac(11, 10) = k9*x6;

jac(11, 11) = - k10 - k11;

jac(11, 20) = x6*x10;

jac(11, 21) = -x11;

jac(11, 22) = -x11;

