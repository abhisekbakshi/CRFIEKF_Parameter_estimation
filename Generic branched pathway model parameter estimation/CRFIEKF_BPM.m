function [M] = CRFIEKF_BPM( time, model, jac, Hk , Q, R, xe0, Pe0, D, d, alpha, T, Ts)

% Detect model size
N=size(Q,2);

% Set initial conditions
xe=zeros(N, (T/Ts));
constrained_xe = xe;

xe(:,1)=xe0;

Pe=cell(1, (T/Ts));
I=eye(N);
Pe{1}=Pe0;

observed = zeros(2,1);

%%% Create Fuzzy inference system
fis_x3_x1 = readfis('x3_x1');
options = evalfisOptions;
options.NumSamplePoints = 200;


fis_x2_x4_x3 = readfis('x2_x4_x3');
options = evalfisOptions;
options.NumSamplePoints = 200;



for ii=2:(T/Ts)
    % **1** Compute linearization around previous a posteriori estimate
    Fk=jac(time(ii-1),xe(:,ii-1));
    
    % **2** A priori estimate of the current state
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol=ode113(model,[time(ii-1) time(ii)], xe(:,ii-1),opts);
    S=size(sol.y);
    xp=sol.y(:,S(2));
    
    % **3** Restrict a priori estimate within range [0,1]
    for j = 1:11
        if (xp(j, 1))<0.1
            xp(j, 1) = 0.1;
        end
    end
    
    for j = 1:11
        if (xp(j, 1))> 0.9
            xp(j, 1) = 0.9;
        end
    end    
         
    
    % **4** Integrates the continuos-time differential Lyapunov equation
    Pp=DifferentialLyapunov113(Fk,Q,[time(ii-1) time(ii)],Pe{ii-1});
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize

    
    % **5** Compute optimal Kalman gain
    Lk=Pp*Hk'*(Hk*Pp*Hk'+R)^(-1);
    
    % **6** Generate measument from Fuzzy inference system
    observed(1,1) = evalfis([xp(3,1)],fis_x3_x1,options);
    observed(2,1) = evalfis([xp(2,1) xp(4,1)],fis_x2_x4_x3,options);
    

    % **7** Incorporate fuzzy measurement to get a posteriori estimate of current prediction 
    xe(:,ii)=xp+Lk*(observed(1:2,1) - Hk*xp);
    Pe{ii}=(I-Lk*Hk)*Pp*(I-Lk*Hk)'+Lk*R*Lk';
    Pe{ii}=(Pe{ii}+Pe{ii}')/2; % Symmetrize
    
    % **8** If condition number > Threshold ==> replace xe with Regularized xe
    reg_mat = pinv(Pe{ii});
    N_k = Hk'*(R^-1)*Hk + Pe{ii};
%     n = norm(N_k);
    con(ii) = norm(N_k) * norm(pinv(N_k));
    
    if (con(ii) > 1000)
        test =  (Hk'*R^(-1)*Hk + Pe{ii} + alpha*reg_mat)^(-1);
        K_trkf = test  * (Hk'*R^(-1));
  
        x_trkf = xp + K_trkf*(observed(1:2,1) - Hk*xp);
  
        P_trkf = (I - K_trkf*Hk)*Pe{ii};
  
        xe(:,ii) = x_trkf(:,1);
        Pe{ii} = P_trkf;
        Pe{ii}=(Pe{ii}+Pe{ii}')/2; 
    end
    
    % **9** If constraints not satisfied ==> replace Regularized xe with constrained regularized xe
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
    
    % **10** Restrict upper bound of constrained regularized xe to 1.
    for j = 1:11
        if (xe(j, ii))> 0.9
            xe(j, ii) = 0.9;
        end
    end    
    

    if(mod(ii, 1000) == 0)
        fprintf("Iteration %d completed\n",ii);
    end
    
end

for i= 1:7
   M(i,1) = mean(xe(i+4,ii-500:ii-1)); 
end


end
