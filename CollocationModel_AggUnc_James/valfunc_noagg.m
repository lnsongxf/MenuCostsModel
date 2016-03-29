function [ v, Phi_pPA ] = valfunc_noagg(flag,cE,s,pPstar,Y,param,glob,options)
%VALFUNC_NOAGG gives value function value given parameters, states

% PREVIOUS: [v,N,Y,I,PI,AC,C,Phi_KpA] = valfunc(flag,c,s,Kp,Y,param,glob,options)
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
    
%-------------------------------------------------

%% Compute flow payoffs if keeping vs. changing price

switch flag
    case 'K'
        PI        = menufun_noagg('PIK',s,[],[],Y,param,glob,options);
        Phi = glob.Phi;
        v          = PI + param.beta*Phi*cE;
    case 'C'
        PI        = menufun_noagg('PIC',s,pPstar,[],Y,param,glob,options);
        Phi_pP  = splibas(glob.pPgrid0,0,glob.spliorder(1),pPstar);
        Phi_pPA     = dprod(glob.Phi_A,Phi_pP);
        v          = PI + param.beta*Phi_pPA*cE;
end

end
