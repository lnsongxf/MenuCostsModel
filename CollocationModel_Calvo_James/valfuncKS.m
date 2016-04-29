function [v,Phi_pPAMY] = valfuncKS(flag,cE,s,pP,param,glob,options)
%VALFUNC gives value function value given parameters, states
%-------------------------------------------------
%   Computes the the value function
%
%   INPUTS
%   - c2         = current collocation coefficient matrix
%   - s         = state space
%   - Kp        = >????
%   - Y         = conjectured value of output, Y
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = value function 
%-------------------------------------------------

switch flag
    case 'K'
        
        % Compute flow payoff
        PI              = menufun('PIK',s,[],param,glob,options);       
        Phi             = funbas(glob.fspace,s);  % glob.Phi;
        v               = PI + param.beta*Phi*cE;
        
    case 'C'
        PI              = menufun('PIC',s,pP,param,glob,options);
        
        % Create basis matrices for continuation value
        Phi_pP          = splibas(glob.pPgrid0,0,glob.spliorder(1),pP);
        Phi_pPAMY       = dprod(glob.Phi_Y, dprod(glob.Phi_M, dprod(glob.Phi_A,Phi_pP)));
        
        % Compute value if changing
        v               = PI + param.beta*Phi_pPAMY*cE;
        
end

end
            