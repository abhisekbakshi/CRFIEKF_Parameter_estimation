function jac = BPM_extended_model_jacobean(t,x_v)

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


% g1 = 0.9;
% g2 = 0.33;
% g3 = 0.05;
% g4 = 0.17;
% g5 = 0.01;
% g6 = 0.22;


g1 = 20;
g2 = 8;
g3 = 3;
g4 = 5;
g5 = 2;
g6 = 6;


% % % g1 = x_v(12);
% % % g2 = x_v(13);
% % % g3 = x_v(14);
% % % g4 = x_v(15);
% % % g5 = x_v(16);
% % % g6 = x_v(17);



x0 = 0.6;


% Preallocate Jacobian
jac=zeros(11,11);
% Set nonzero elements

jac(1, 1) = - f21*g2*x1^(f21 - 1) - f51*g5*x1^(f51 - 1);

jac(1, 3) = -(f13*g1*x0)/x3^(f13 + 1);

jac(1, 5) = -(g1*x0*log(x3))/x3^f13;

jac(1, 6) = -g2*x1^f21*log(x1);

jac(1, 10) = -g5*x1^f51*log(x1);

jac(2, 1) = f21*g2*x1^(f21 - 1);

jac(2, 2) = -f32*g3*x2^(f32 - 1);

jac(2, 6) = g2*x1^f21*log(x1);

jac(2, 7) = -g3*x2^f32*log(x2);

jac(3, 2) = f32*g3*x2^(f32 - 1);

jac(3, 3) = -f43*g4*x3^(f43 - 1)*x4^f44;

jac(3, 4) = -f44*g4*x3^f43*x4^(f44 - 1);

jac(3, 7) = g3*x2^f32*log(x2);

jac(3, 8) = -g4*x3^f43*x4^f44*log(x3);

jac(3, 9) = -g4*x3^f43*x4^f44*log(x4);

jac(4, 1) = f51*g5*x1^(f51 - 1);

jac(4, 4) = -f64*g6*x4^(f64 - 1);

jac(4, 10) = g5*x1^f51*log(x1);

jac(4, 11) = -g6*x4^f64*log(x4);




% % % % Preallocate Jacobian
% % % jac=zeros(17,17);
% % % % Set nonzero elements
% % % 
% % % jac(1, 1) = - f21*g2*x1^(f21 - 1) - f51*g5*x1^(f51 - 1);
% % % 
% % % jac(1, 3) = -(f13*g1*x0)/x3^(f13 + 1);
% % % 
% % % jac(1, 5) = -(g1*x0*log(x3))/x3^f13;
% % % 
% % % jac(1, 6) = -g2*x1^f21*log(x1);
% % % 
% % % jac(1, 10) = -g5*x1^f51*log(x1);
% % % 
% % % jac(1, 12) = x0/x3^f13;
% % % 
% % % jac(1, 13) = -x1^f21;
% % % 
% % % jac(1, 16) = -x1^f51;
% % % 
% % % jac(2, 1) = f21*g2*x1^(f21 - 1);
% % % 
% % % jac(2, 2) = -f32*g3*x2^(f32 - 1);
% % % 
% % % jac(2, 6) = g2*x1^f21*log(x1);
% % % 
% % % jac(2, 7) = -g3*x2^f32*log(x2);
% % % 
% % % jac(2, 13) = x1^f21;
% % % 
% % % jac(2, 14) = -x2^f32;
% % % 
% % % jac(3, 2) = f32*g3*x2^(f32 - 1);
% % % 
% % % jac(3, 3) = -f43*g4*x3^(f43 - 1)*x4^f44;
% % % 
% % % jac(3, 4) = -f44*g4*x3^f43*x4^(f44 - 1);
% % % 
% % % jac(3, 7) = g3*x2^f32*log(x2);
% % % 
% % % jac(3, 8) = -g4*x3^f43*x4^f44*log(x3);
% % % 
% % % jac(3, 9) = -g4*x3^f43*x4^f44*log(x4);
% % % 
% % % jac(3, 14) = x2^f32;
% % % 
% % % jac(3, 15) = -x3^f43*x4^f44;
% % % 
% % % jac(4, 1) = f51*g5*x1^(f51 - 1);
% % % 
% % % jac(4, 4) = -f64*g6*x4^(f64 - 1);
% % % 
% % % jac(4, 10) = g5*x1^f51*log(x1);
% % % 
% % % jac(4, 11) = -g6*x4^f64*log(x4);
% % % 
% % % jac(4, 16) = x1^f51;
% % % 
% % % jac(4, 17) = -x4^f64;


end






