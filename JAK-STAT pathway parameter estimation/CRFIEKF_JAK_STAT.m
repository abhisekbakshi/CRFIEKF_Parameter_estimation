function [M] = CRFIEKF_JAK_STAT( time, model, jac, Hk , Q, R, xe0, Pe0, D, d, alpha)


% Detect model size
N=size(Q,2);

% Set initial conditions
xe=zeros(N, 16);
constrained_xe = xe;

xe(:,1)=xe0;

Pe=cell(1,16);
I=eye(N);
Pe{1}=Pe0;


%%% Create Fuzzy inference system
fis = readfis('measurement');
options = evalfisOptions;
options.NumSamplePoints = 200;

for ii=2:5000
    % **1** Compute linearization around previous a posteriori estimate
    Fk=jac(time(ii-1),xe(:,ii-1));
    
    % **2** A priori estimate of the current state
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol=ode113(model,[time(ii-1) time(ii)], xe(:,ii-1),opts);
    S=size(sol.y);
    xp=sol.y(:,S(2));
    
    y1 = xp(2) + 2*xp(3);
    y2(ii-1) = xp(1) + xp(2) + 2*xp(3);


    % **3** Integrates the continuos-time differential Lyapunov equation
    Pp=DifferentialLyapunov113(Fk,Q,[time(ii-1) time(ii)],Pe{ii-1});
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize

    
    % **4** Compute optimal Kalman gain
    Lk=Pp*Hk'*(Hk*Pp*Hk'+R)^(-1);
    
    % **5** Generate measument from Fuzzy inference system
    observed = evalfis(y1,fis,options);
    ob(ii-1) = observed;

    % **6** Incorporate fuzzy measurement to get a posteriori estimate of current prediction 
    xe(:,ii)=xp+Lk*(observed - y2(ii-1));
    Pe{ii}=(I-Lk*Hk)*Pp*(I-Lk*Hk)'+Lk*R*Lk';
    Pe{ii}=(Pe{ii}+Pe{ii}')/2; % Symmetrize
    
    % **7** If condition number > Threshold ==> replace xe with Regularized xe
    reg_mat = pinv(Pe{ii});
    N_k = Hk'*(R^-1)*Hk + Pe{ii};
%     n = norm(N_k);
    con(ii) = norm(N_k) * norm(pinv(N_k));
    
    if (con(ii) > 2446.47)
        test =  (Hk'*R^(-1)*Hk + Pe{ii} + alpha*reg_mat)^(-1);
        K_trkf = test  * (Hk'*R^(-1));
  
        x_trkf = xp + K_trkf*(observed - Hk*y2(ii-1));
  
        P_trkf = (I - K_trkf*Hk)*Pe{ii};
  
        xe(:,ii) = x_trkf(:,1);
        Pe{ii} = P_trkf;
        Pe{ii}=(Pe{ii}+Pe{ii}')/2; 
    end
    
    % **8** If constraints not satisfied ==> replace Regularized xe with constrained regularized xe
    if ~all(D*xe(:,ii)<=d)
        if all(eig(Pe{ii})>zeros(size(Pe{ii},1),1))
            W=inv(Pe{ii});
            W=(W+W')/2;
        else
            disp('Warning. Non positive-definite Pe.')
            W=eye(size(Pe{ii},1));
            break;
        end
        % Solve the quadratic programming problem
        oo = optimset('Algorithm','interior-point-convex','Display','off');
        z = quadprog(2*W, zeros(size(xe(:,ii-1))), D, d-D*xe(:,ii),[],[],[],[],[],oo);
        
        % Replace unconstrained estimate with constrained estimate
        constrained_xe(:,ii) = xe(:,ii) + z;
        xe(:,ii) = constrained_xe(:,ii);
    else
        xe(:,ii) = xe(:,ii);
    end
    
    if(mod(ii, 1000) == 0)
        fprintf("Iteration %d completed\n",ii);
    end
     
    
end

%%% Estimating result
for i= 1:4
   M(i,1) = mean(xe(i+4,ii-200:ii-1)); 
end


end
