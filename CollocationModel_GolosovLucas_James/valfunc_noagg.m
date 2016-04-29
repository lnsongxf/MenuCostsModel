function [ v, Phi ] = valfunc_noagg(flag,cE,s,pPstar,Y,param,glob,options)
%VALFUNC_NOAGG gives value function value given parameters, states
%-------------------------------------------------
%   Computes the the value function
%
%   INPUTS
%   - c2         = current collocation coefficient matrix
%   - s         = state space
%   - pPstar    = Optimal real price choice of the firm
%   - Y         = conjectured value of output, Y
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = value function 
%   - Phi_pPV   = basis matrix given new price choice, pPstar.
    
%-------------------------------------------------

%% Compute flow payoffs if keeping vs. changing price

%
% NOTE: Should we be updating tomorrow's price inside the expectation term?
% I think we have been missing this previously...
%

switch flag
    case 'K'
        PI          = menufun('PIK',s,[],[],Y,param,glob,options);
        s_prime     = [s(:,1)*(1/exp(glob.piw)), s(:,2)];
        Phi         = funbas(glob.fspace,s);    % s_prime);
        v           = PI + param.beta*Phi*cE;
    case 'C'
        PI          = menufun('PIC',s,pPstar,[],Y,param,glob,options);
        Phi_pP      = splibas(glob.pPgrid0,0,glob.spliorder(1),pPstar); %*1/exp(glob.piw)); 
        Phi_V       = splibas(glob.vgrid0,0,glob.spliorder(2),s(:,2));
        Phi         = dprod(Phi_V,Phi_pP);
        v           = PI + param.beta*Phi*cE;
end

end
