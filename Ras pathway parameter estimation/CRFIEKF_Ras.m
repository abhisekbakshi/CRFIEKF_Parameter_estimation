function [M] = CRFIEKF_Ras(time, model, jac, Hk , Q, R, xe0, Pe0, alpha)

% Detect model size
N=size(Q,2);

% Set initial conditions
xe=zeros(N, 11);

xe(:,1)=xe0;

Pe=cell(1,11);
I=eye(N);
Pe{1}=Pe0;


%%% Create Fuzzy inference system
fis = readfis('measurement');
options = evalfisOptions;
options.NumSamplePoints = 200;


for ii=2:11
    % **1** Compute linearization around previous a posteriori estimate
    Fk=jac(time(ii-1),xe(:,ii-1));
    
    % **2** A priori estimate of the current state
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol=ode113(model,[time(ii-1) time(ii)], xe(:,ii-1),opts);
    S=size(sol.y);
    xp=sol.y(:,S(2));
    
    % **3** Integrates the continuos-time differential Lyapunov equation
    Pp=DifferentialLyapunov113(Fk,Q,[time(ii-1) time(ii)],Pe{ii-1});
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize
    
    % **4** Compute optimal Kalman gain
    Lk=Pp*Hk'*(Hk*Pp*Hk'+R)^(-1);
    
    % **5** Generate measument from Fuzzy inference system
    observed(:,1) = evalfis([xp(7), xp(10)],fis,options);
    
    % **6** Incorporate fuzzy measurement to get a posteriori estimate of current prediction 
    xe(:,ii)=xp+Lk*(observed(:,1) - Hk*xp);
    Pe{ii}=(I-Lk*Hk)*Pp*(I-Lk*Hk)'+Lk*R*Lk';
    Pe{ii}=(Pe{ii}+Pe{ii}')/2; % Symmetrize
    
    % **7** If condition number > Threshold ==> replace xe with Regularized xe
    reg_mat = pinv(Pe{ii});
    N_k = Hk'*(R^-1)*Hk + Pe{ii};
    con(ii) = norm(N_k) * norm(pinv(N_k));
    
    if (con(ii) > 2.7179e+06)
        test =  (Hk'*R^(-1)*Hk + Pe{ii} + alpha*reg_mat)^(-1);
        K_trkf = test  * (Hk'*R^(-1));
  
        x_trkf = xp + K_trkf*(observed(:,1) - Hk*xp);
  
        P_trkf = (I - K_trkf*Hk)*Pe{ii};
  
        xe(:,ii) = x_trkf(:,1);
        Pe{ii} = P_trkf;
        Pe{ii}=(Pe{ii}+Pe{ii}')/2; 
    end
   
end

%%% Result
M = xe(12:22, 10);

end

