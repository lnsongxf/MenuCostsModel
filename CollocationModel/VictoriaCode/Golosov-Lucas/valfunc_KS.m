function [v,Phi_XpNu] = valfunc_KS(c,s,xp,param,glob,options)

    %__________________________________________________________________________
    % Compute flow profits if change
    x           = s(:,1);
    nu          = s(:,2);
    cbar        = s(:,3);
    F           = menufun_KS('change',[x,nu,cbar],xp,param,glob,options);
    %__________________________________________________________________________
    % Create basis matrices for continuation value
    Phi_Xp      = splibas(glob.xgrid0,0,glob.spliorder(1),xp);
    Phi_XpNu    = dprod(glob.Phi_c,dprod(glob.Phi_nu,Phi_Xp));    
    %__________________________________________________________________________
    % Compute value
    ce         = c(2*end/3+1:end);
    v          = F + param.beta*Phi_XpNu*ce;       

end