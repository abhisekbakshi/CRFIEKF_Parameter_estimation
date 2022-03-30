clear all;
close all;
clc;


x1 = 0.7;
x2 = 0.65;
x3 = 0.52;
x4 = 0.75;

x_v = zeros(1,4);
x_v(1,:) = [x1 x2 x3 x4];


T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';

for ii=2:5000
    [tsim, xsim] = ode23tb('BPM_model', Tspan, x_v(ii-1,:));
    x_v(ii, :) = real(xsim(end, :));
    
end

