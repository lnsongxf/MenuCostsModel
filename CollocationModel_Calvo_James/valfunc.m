function [ v, Phi_pPAMY ] = valfunc(flag,cE,s,pPstar,Y,param,glob,options)
%VALFUNC gives value function value given parameters, states

% PREVIOUS: [v,N,Y,I,PI,AC,C,Phi_KpA] = valfunc(flag,c,s,Kp,Y,param,glob,options)
%-------------------------------------------------
%   Computes the the value function
%
%   INPUTS
%   - c2        = current collocation coefficient matrix
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

% pP          = s(:,1);
% A           = s(:,2);
switch flag
    case 'K'
        PI          = menufun('PIK',s,[],Y,param,glob,options); 
        Phi         = glob.Phi;
        v           = PI + param.beta*Phi*cE;  
        
    case 'C'
        PI        = menufun('PIC',s,pPstar,Y,param,glob,options); 
        Phi_pP  = splibas(glob.pPgrid0,0,glob.spliorder(1),pPstar);
        Phi_pPAMY     = dprod(glob.Phi_Y, dprod(glob.Phi_M, dprod(glob.Phi_A,Phi_pP)));    
        v          = PI + param.beta*Phi_pPAMY*cE;       

end
   
%__________________________________________________________________________
%% If requested, provide other policies
% if nargout>1    
%     N       = menufun('labor',[pP,A],Kp,P,param,glob,options);
%     Y       = menufun('output',[pP,A],Kp,P,param,glob,options);
%     I       = menufun('investment',[pP,A],Kp,P,param,glob,options);
%     AC      = menufun('costs',[pP,A],Kp,P,param,glob,options);
%     C       = Y - I - AC;
% end   




