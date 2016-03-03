function [v,jac] = solve_valfunc_menucost(ce,s,Y,param,glob,options)
%__________________________________________________________________________
% First, use golden search to solve for price if the firm changes (Pc): 
obj                     = @(pc)valfuncchange(ce,s,pc,Y,param,glob,options);
Pc                      = goldenx(obj,ones(size(s,1),1)*0,ones(size(s,1),1)*5);
% vc is the value if the firm changes:
[vc,Phi_PcA]          = valfuncchange(ce,s,Pc,Y,param,glob,options);

% Next, compute the value if the firm keep its price
Pikeep = menufun_menucosts('keep',s,s(:,1),Y,param,glob,options);
vk     = Pikeep + param.beta*funbas(glob.fspace,s)*ce;

% Find the maximum of vc and vk and define I(s)
maxval = max(vk,vc);
Is = double((vc>vk));

% Find the final RHS of value function
% to account for other grid sizes, need to find length of p grid
pgrid   = s(s(:,2)==s(1,2),1); 
plength = size(pgrid,1);
vf = kron(glob.A,speye(plength))*maxval;

% Find a policy function for prices
Pp = Pc.*Is + s(:,1).*(1-Is);

%__________________________________________________________________________
% Compute jacobian if requested
if (nargout==2)
    jac = funbas(glob.fspace,s) - param.beta*kron(glob.A,speye(plength))...
        *(dprod((1-Is),funbas(glob.fspace,s)) + dprod(Is,Phi_PcA));
end
%__________________________________________________________________________
% Packup output
v.vf    = vf;
v.vc    = vc;
v.vk    = vk;
v.Pc    = Pc;
v.Pp    = Pp;
v.Is    = Is;

end