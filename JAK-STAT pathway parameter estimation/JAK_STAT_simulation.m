function [x_v, ph_stat5_observed_scaled, total_stat5_observed_scaled] = JAK_STAT_simulation(estimated_parameter)


% % clear all;
% % close all;
% % clc;

global k1 k2 k3 k4;

k1 = estimated_parameter(1);
k2 = estimated_parameter(2);
k3 = estimated_parameter(3);
k4 = estimated_parameter(4);

u = [0,9,14,44,59,50,46,32,35,19,19,17,3,2,1,1];
time = [0,2,4,6,8,10,12,14,16,18,20,25,30,40,50,60];
x_v = zeros(16,5);
x_v(1,:) = [0.5 0 0 0 0];

ph_stat5_observed = zeros(1,16);
total_stat5_observed = zeros(1,16);

for i=2:16
    [tsim, xsim] = ode23tb('JAK_STAT_model', [time(i-1) time(i)], x_v(i-1,:));
    x_v(i, 1:4) = real(xsim(end, 1:4));
    x_v(i, 5) = u(i);
    
    ph_stat5_observed(1,i) = x_v(i,2) + 2*x_v(i,3);
    total_stat5_observed(1,i) = x_v(i,1) + x_v(i,2) + 2*x_v(i,3);
    
end

ph_stat5_observed_scaled = rescale(ph_stat5_observed);
total_stat5_observed_scaled = rescale(total_stat5_observed);

end

