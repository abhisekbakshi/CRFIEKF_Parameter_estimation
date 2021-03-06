%%% DifferentialLyapunov

% Integrates the continuos-time differential Lyapunov equation
%
%       Pdot=F*P+P*F'+Q

function P=DifferentialLyapunov113(F,Q,tspan,P0)

p0=reshape(P0,length(F)^2,1);
opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
[tsim,psim]=ode113(@(t,y)dreof(t,y,F,Q),tspan,p0,opts);

P=cell(length(tsim),1);

for ii=1:length(tsim)
    P{ii}=reshape(psim(ii,:),size(F));
end


%%% Ode Function for Differential Lyapunov

function yp=dreof(t,y,F,Q)

Y=reshape(y,size(F));
Y=(Y+transpose(Y))./2; % Symmetrize
YP=F*Y+Y*transpose(F)+Q;

yp=reshape(YP,length(F)^2,1);
