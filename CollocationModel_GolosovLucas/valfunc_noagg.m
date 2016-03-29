function [ v, Phi_pPV ] = valfunc_noagg(flag,cE,s,pPstar,Y,param,glob,options)
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

switch flag
    case 'K'
        PI         = menufun('PIK',s,[],[],Y,param,glob,options);
        Phi        = glob.Phi;
        v          = PI + param.beta*Phi*cE;
    case 'C'
        PI        = menufun('PIC',s,pPstar,[],Y,param,glob,options);
        Phi_pP    = splibas(glob.pPgrid0,0,glob.spliorder(1),pPstar);
        Phi_pPV   = dprod(glob.Phi_V,Phi_pP);
        v         = PI + param.beta*Phi_pPV*cE;
end

end
