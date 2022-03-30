
%%%
% Function version of the extended Kalman filter for parameter estimation

% Inputs
% time  ==> time vector
% model ==> function handle to the model used for estimation
% M     ==> measurements vector

% Outputs are
% xe    ==> estimated state
% Pe    ==> percentual estimation error

function [xe,Pe,R]=parameter_estimation_filter_HEKF(time,model,jac,M,Hk,Q,R,xe0,Pe0,D,d)

%%% HEKF run

% Detect model size
N=size(Q,2);

% **0** Set initial conditions
xe=zeros(N, size(M,2));
% disp(size(xe,1));
xe(:,1)=xe0;



Pe=cell(1,size(M,2));
I=eye(N);
Pe{1}=Pe0;

bar=waitbar(0,'Running HEKF - Time = 0');

for ii=2:length(M)
    % **1** Compute linearization around previous a posteriori estimate
    Fk=jac(time(ii-1),xe(:,ii-1));
    
    % **2** Advance time
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol=ode113(model,[time(ii-1) time(ii)], xe(:,ii-1),opts);
    S=size(sol.y);
    xp=sol.y(:,S(2));
    
    for j = 1:71
        if (xp(j, 1))<0.1
            xp(j, 1) = 0.1;
        end
    end
   
    for j = 1:71
        if (xp(j, 1))>1.0
            xp(j, 1) = 1.0;
        end   
    end
    
    
    % disp(xp)

    Pp=DifferentialLyapunov113(Fk,Q,[time(ii-1) time(ii)],Pe{ii-1});
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize
%      disp(Pp)
    
    % **3** Compute optimal gain (try to avoid inverse...)
    Lk=Pp*Hk'*(Hk*Pp*Hk'+R)^(-1);
    % Lk=linsolve(Hk*Pp*Hk'+R, Hk*Pp);
    % Lk=Lk';
    % disp(Lk)
%     error(:,ii) = M(:,ii)-Hk*xp;
    % **4** Incorporate currente measurement
    xe(:,ii)=xp+Lk*(M(:,ii)-Hk*xp);
    Pe{ii}=(I-Lk*Hk)*Pp*(I-Lk*Hk)'+Lk*R*Lk';
    Pe{ii}=(Pe{ii}+Pe{ii}')/2; % Symmetrize
    

    
    % disp(xe(:,ii));
    
    % **5** If constraints not satisfied ==> replace xe with constrained estimate
    if ~all(D*xe(:,ii)<=d)
%         disp('Warning. Constraints not satified. Recomputing constrained estimate.')
        % Check invertibility and positive definiteness of Pe{ii}
        if all(eig(Pe{ii})>zeros(size(Pe{ii},1),1))
            W=inv(Pe{ii});
            W=(W+W')/2;
        else
            disp('Warning. Non positive-definite Pe.')
            W=eye(size(Pe{ii},1));
        end
        % Solve the quadratic programming problem
        oo=optimset('Display', 'off','MaxFunEvals',1e4, 'MaxIter',1e2, 'LargeScale', 'off');
        z = quadprog(2*W, zeros(size(xe(:,ii))), D, d-D*xe(:,ii),[],[],zeros(71,1),ones(71,1),[],oo);
        
        % Replace unconstrained estimate with constrained estimate
        xe(:,ii) = xe(:,ii) + z;
    end
    
    waitbar(ii/length(M),bar,['Running HEKF - Time = ' num2str(time(ii))]);
end

%%% 

close(bar);

end