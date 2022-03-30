function jac = JAK_STAT_model_jacobian(t,y,u)


x1 = y(1);
x2 = y(2);
x3 = y(3);
x4 = y(4);

k1 = y(5);
k2 = y(6);
k3 = y(7);
k4 = y(8);

u = (2*k4*x3)/(k1*x1);



% Preallocate Jacobian
jac=zeros(8,8);
% Set nonzero elements

jac(1, 1) = -k1*u;

jac(1, 3) = 2*k4;

jac(1, 5) = -u*x1;

jac(1, 8) = 2*x3;

jac(2, 1) = k1*u;

jac(2, 2) = -4*k2*x2;

jac(2, 5) = u*x1;

jac(2, 6) = -2*x2^2;

jac(3, 2) = 2*k2*x2;

jac(3, 3) = -k3;

jac(3, 6) = x2^2;

jac(3, 7) = -x3;

jac(4, 3) = k3 - k4;

jac(4, 7) = x3;

jac(4, 8) = -x3;

end






