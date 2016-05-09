function [vc,Phi_PpA] = valfuncchange(ce,s,pp,Y,param,glob,options)
%__________________________________________________________________________
% Compute flow profits if change
p           = s(:,1);
a           = s(:,2);
F           = menufun_menucosts('change',[p,a],pp,Y,param,glob,options);
%__________________________________________________________________________
% Create basis matrices for continuation value
Phi_Pp      = splibas(glob.pgrid0,0,glob.spliorder(1),pp*1/exp(param.mu));
Phi_PpA     = dprod(glob.Phi_A,Phi_Pp);    
%__________________________________________________________________________
% Compute value
vc          = F + param.beta*Phi_PpA*ce;       
%__________________________________________________________________________
% % If requested, provide other policies
% if nargout>1
%     N       = menufun('labor',[K,Z],Kp,P,param,glob,options);
%     Y       = menufun('output',[K,Z],Kp,P,param,glob,options);
%     I       = menufun('investment',[K,Z],Kp,P,param,glob,options);
%     AC      = menufun('costs',[K,Z],Kp,P,param,glob,options);
%     C       = Y - I - AC;
% end   
end