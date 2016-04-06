function [v,Phi_XpNu] = valfunc(c,s,xp,cbar,param,glob,options)

    %__________________________________________________________________________
    % Compute flow profits if change
    x           = s(:,1);
    nu          = s(:,2);
    F           = menufun('flow',[x,nu],xp,cbar,param,glob,options);
    %__________________________________________________________________________
    % Create basis matrices for continuation value
    Phi_Xp      = splibas(glob.xgrid0,0,glob.spliorder(1),xp);
    Phi_XpNu    = dprod(glob.Phi_Nu,Phi_Xp);    
    %__________________________________________________________________________
    % Compute value
    v          = F + param.beta*Phi_XpNu*ce;       

end