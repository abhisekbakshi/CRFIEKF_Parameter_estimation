function [x_v, ref] = BPM_model_simulation(k)

global f13 f21 f32 f43 f44 f51 f64;


f13 = k(1);
f21 = k(2);
f32 = k(3);
f43 = k(4);
f44 = k(5);
f51 = k(6);
f64 = k(7);

x1 = 1.5;
x2 = 2.9;
x3 = 1.1;
x4 = 0.47;



x_v = zeros(1,4);
x_v(1:4) = [x1 x2 x3 x4];

T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';

for ii = 1:(T/Ts)
                                                            
 [tsim, xsim] = ode23tb('BPM_model', Tspan, x_v(ii,:));
 x_v(ii+1, :) = real(xsim(end, :));

end
disp('Generic branch pathway model simulation using estimated parameter completed...');



ref = zeros(1,4);
ref(1:4) = [x1 x2 x3 x4];

T  = 10;
Ts = 0.002;
Tspan = [0 Ts];
t = [0:Ts:T]';

for ii = 1:(T/Ts)
                                                            
 [tsim, xsim] = ode23tb('BPM_model_original', Tspan, ref(ii,:));
 ref(ii+1, :) = real(xsim(end, :));

end

disp('Generic branch pathway model simulation using previously estimated parameter completed...');



end