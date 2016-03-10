function [v1,N,Y,I,F,AC,C,Phi_kKzA] = valfuncKS(c2,s,kp,P,param,glob,options)
%__________________________________________________________________________
% Compute flow payoff
k           = s(:,1);
z           = s(:,3); 
A           = s(:,4); 
F           = menufun('F',[k,z.*A],kp,P,param,glob,options);   
%__________________________________________________________________________
% Create basis matrices for continuation value
Phi_k       = splibas(glob.kgrid0,0,glob.spliorder(1),kp);
Phi_kKzA    = dprod(glob.Phi_KzA,Phi_k); 
%__________________________________________________________________________
% Compute value
v1          = F + glob.beta*Phi_kKzA*c2;               
%__________________________________________________________________________
% If requested, provide other policies
if nargout>1
    N       = menufun('labor',[k,z.*A],kp,P,param,glob,options);
    Y       = menufun('output',[k,z.*A],kp,P,param,glob,options);
    I       = menufun('investment',k,kp,P,param,glob,options);
    AC      = menufun('costs',k,kp,P,param,glob,options);
    C       = Y - I - AC;
end 