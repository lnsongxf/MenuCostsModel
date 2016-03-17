function [v,jac] = solve_valfuncKS(c,s,param,glob,options,xxx)

ck = c(1:end/3);  
cc = c(end/3+1:2*end/3);
ce = c(2*end/3+1:end);

%__________________________________________________________________________
% First, use golden search to solve for price if the firm changes (Pc): 
obj                     = @(pc)valfuncchange_KS(ce,s,pc,param,glob,options);
Pc                      = goldenx(obj,ones(size(s,1),1)*min(glob.pgrid),ones(size(s,1),1)*max(glob.pgrid));
% vc is the value if the firm changes:
[vc,Phi_PpA]          = valfuncchange_KS(ce,s,Pc,param,glob,options);

% Next, compute the value if the firm keep its price
Pikeep = menufun_menucosts_KS('keep',s,s(:,1),param,glob,options);
vk     = Pikeep + param.beta*funbas(glob.fspace,s)*ce;

% Compute I(s)
Is = double((vc>vk));
%Is = double((glob.Phi_stilde*cc>=glob.Phi_stilde*ck));

% policy function
Pp = Pc.*Is + s(:,1).*(1-Is);
%Pp = kron(Pc,ones(glob.Npi,1)).*Is + kron(s(:,1),ones(glob.Npi,1)).*(1-Is);
Is = kron(Is,ones(glob.Npi,1));

% RHS of expected value equation
%ve = glob.exp_matrix*(max(glob.Phi_stilde*ck,glob.Phi_stilde*cc));
ve  = [];
if (nargin<=5)  % When simulating etc, don't need this, nargin is useful
    ve = glob.exp_matrix*(dprod((1-Is),glob.Phi_stilde)*ck + dprod(Is,glob.Phi_stilde)*cc);     
end

%__________________________________________________________________________
% Compute jacobian if requested
if (nargout==2)
    jac         = [glob.Phi, sparse(glob.Ns,glob.Ns), -param.beta*glob.Phi; ...
                   sparse(glob.Ns,glob.Ns), glob.Phi, -param.beta*Phi_PpA;
                   -glob.exp_matrix*(dprod((1-Is),glob.Phi_stilde)), ...
                   -glob.exp_matrix*(dprod(Is,glob.Phi_stilde)), glob.Phi];
                   
end
%__________________________________________________________________________
% Packup output
v.vc    = vc;
v.vk    = vk;
v.ve    = ve;
v.Pc    = Pc;
v.Pp    = Pp;
v.Is    = Is;

end