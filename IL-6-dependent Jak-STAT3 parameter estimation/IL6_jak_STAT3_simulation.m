function [x_v, reference] = IL6_jak_STAT3_simulation(k)


% clc
% clear
% close all

warning('off','all');

global p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13;

p1	=	k(1);
p2	=	k(2);
p3	=	k(3);
p4	=	k(4);
p5	=	k(5);
p6	=	k(6);
p7	=	k(7);
p8	=	k(8);
p9	=	k(9);
p10	=	k(10);
p11	=	k(11);
p12	=	k(12);
p13	=	k(13);


%%%%% METABOLITES INITIAL CONCENTRATION in nM %%%%%

x1	=	0;     %	IL-6IL-6R
x2	=	16.8;  %	gp130
x3	=	0;     %	(p)Rcomplex
x4	=	0;     %	(p)STAT3
x5	=	0;     %	SOCS3 mRNA
x6	=	0;     %	SOCS3 protein


global x8;
global x7;
global x9;

x8 = (16.8 - x2 - 2*x3)/2;      %  Rcomplex
x7 = 2.1 - x1 - 2*x3 - 2*x8;    %  IL-6R
x9 = 83 - x4;                   %  STAT3

reference = zeros(1,6);
all = zeros(1,9);

reference(1:6) = [x1 x2 x3 x4 x5 x6];
all(1:9) = [x1 x2 x3 x4 x5 x6 x7 x8 x9];


T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';
tic

for ii = 1:(T/Ts)
                                    
    [tsim, xsim] = ode23tb('IL6_jak_STAT3_model_original', Tspan, reference(ii,:));
    reference(ii+1, :) = real(xsim(end, :));    
    all(ii+1, :) = [reference(ii+1, :) x7 x8 x9];

end

disp('Simulation of IL6 dependent JAK/STAT pathway using original model parameter completed...');    





x_v = zeros(1,6);
all = zeros(1,9);

x_v(1:6) = [x1 x2 x3 x4 x5 x6];
all(1:9) = [x1 x2 x3 x4 x5 x6 x7 x8 x9];


% uu = zeros(1, 27); 

T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';
tic

for ii = 1:(T/Ts)
                                    
    [tsim, xsim] = ode23tb('IL6_jak_STAT3_model', Tspan, x_v(ii,:));
    x_v(ii+1, :) = real(xsim(end, :));
    all(ii+1, :) = [x_v(ii+1, :) x7 x8 x9];

end

disp('Simulation of IL6 dependent JAK/STAT pathway using estimated parameter completed...');



end



