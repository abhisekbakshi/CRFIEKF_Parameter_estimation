function [x_v] = Ras_model_simulation(k)

global k1 k2 k3 k4 k5 k6 k7 k8 k9 k10  k11;

k1=	k(1);
k2=	k(2);
k3=	k(3);
k4=	k(4);
k5=	k(5);
k6=	k(6);
k7=	k(7);
k8=	k(8);
k9=	k(9);
k10= k(10);
k11= k(11);


time = [0 1 2 3 4 5 6 7 8 9 10];
x_v = zeros(11,11);
x_v(1,:) = [120	48	1.5	4.8	12.4	4.8	87.8	3.7	240.2	160.5	2.7];

for i=2:11
    [tsim, xsim] = ode23tb('Ras_model_for_validation', [time(i-1) time(i)], x_v(i-1,:));
    x_v(i, 1:11) = real(xsim(end, 1:11));   
end


end

