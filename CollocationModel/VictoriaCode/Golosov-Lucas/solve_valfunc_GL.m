function [v,jac] = solve_valfunc_GL(c,s,cbar,param,glob,options,xxx)

    % Extract coefficients
    ck    = c(1:end/3);
    cc    = c(end/3+1:2*end/3);
    ce    = c(2*end/3+1:end);

    %__________________________________________________________________________
    % First, use golden search to solve for price if the firm changes (Xc): 
    obj                     = @(xc)valfunc(c,s,xc,cbar,param,glob,options);
    Xc                      = goldenx(obj,ones(size(s,1),1)*min(glob.xgrid),ones(size(s,1),1)*max(glob.xgrid));
    % vc is the value if the firm changes:
    [vc,Phi_XpNu]           = valfunc(c,s,Xc,cbar,param,glob,options);
    vc                      = vc - param.k;

    % Next, compute the value if the firm keep its price
    Pikeep = menufun('keep',s,s(:,1),cbar,param,glob,options);
    vk     = Pikeep + param.beta*funbas(glob.fspace,s)*ce;

    % Find the maximum of vc and vk and define I(s)
    maxval = max(vk,vc);
    Is = double((vc>vk));

    if (nargin<=6)
        % Find the RHS of expected value function
        % to account for other grid sizes, need to find length of x grid
        xgrid   = s(s(:,2)==s(1,2),1); 
        xlength = size(xgrid,1);
        ve = kron(glob.Nu,speye(xlength))*maxval;
        v.ve    = ve;
    end

    % Find a policy function for prices
    Xp = Xc.*Is + s(:,1).*(1-Is);

    %__________________________________________________________________________
    % Compute jacobian if requested
    if (nargout==2)
        jac     = [funbas(glob.fspace,s), sparse(glob.Ns,glob.Ns), -param.beta*funbas(glob.fspace,s); ...
             	sparse(glob.Ns,glob.Ns), funbas(glob.fspace,s), -param.beta*Phi_XpNu;
                -kron(glob.Nu,speye(glob.Nx))*(dprod((1-Is),glob.Phi)), ...
                -kron(glob.Nu,speye(glob.Nx))*(dprod(Is,glob.Phi)), funbas(glob.fspace,s)];
    end
    %__________________________________________________________________________

    % Packup output
    v.vc    = vc;
    v.vk    = vk;
    v.Xc    = Xc;
    v.Xp    = Xp;
    v.Is    = Is;

end