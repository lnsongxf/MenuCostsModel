function out = menufun(flag,s,pPstar,Y,param,glob,options)
%MENUFUN 
%-------------------------------------------------
%
%   INPUTS
%   - flag      = 
%   - s         = state space
%   - x         = 
%   - Y         = conjectured value of output, Y
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = 
%   - jac       = Jacobian of the value functions  
%-------------------------------------------------

%%
% Globals
    % None
% Parameters
delta   = param.delta;% relative weight of labor-to-consumption in utility
sigma   = param.sigma;% risk aversion coefficient
phielas = param.phielas;% inveser labour supply elasticity
theta   = param.theta;% elasticity of substitution
alpha   = param.alpha;% returns to labour
Phicost = param.Phicost;% menu cost in labour units

%% Cases
switch flag
    case 'bounds'
        pPpmin   = ones(size(s,1),1)*min(glob.pPgrid); 
        pPpmax   = ones(size(s,1),1)*max(glob.pPgrid);
        out     = [pPpmin,pPpmax];
    case 'PIK'
        pP      = s(:,1);
        A       = s(:,2);
        out = Y.*(pP).^(1-theta) - ...
            delta*Y.^sigma.*( (Y./A).^(1/alpha).*(pP).^(-theta/alpha) ).^(1+phielas);
    case 'PIC'
        A       = s(:,2);
        MenuCost = MCfun();
        out = Y.*(pPstar).^(1-theta) - ...
            delta*Y.^sigma.*( (Y./A).^(1/alpha).*(pPstar).^(-theta/alpha) + MenuCost ).^(1+phielas);      
%     case 'output'
%         K       = s(:,1);
%         Z       = s(:,2);
%         N       = Nfun(Z,K);
%         out     = Yfun(Z,N,K); 
%     case 'labor'
%         K       = s(:,1);
%         Z       = s(:,2);
%         out     = Nfun(Z,K);
%     case 'investment'
%         K       = s(:,1);
%         Kp      = x;
%         out     = Ifun(K,Kp);
%     case 'costs'
%         K       = s(:,1);
%         Kp      = x;
%         out     = ACfun(K,Kp);  
end


%% Equilibrium variables and their dependents
% w           = psi./Y;


%__________________________________________________________________________
% NESTED FUNCTIONS

% 1. Menu costs
function MenuCost = MCfun()
    switch options.MC
        case 'Y'
            MenuCost = Phicost*ones(size(s,1),1); 
        case 'N'
            MenuCost = zeros(size(s,1),1);
    end
end

% 2. Labor demand from FOC(n')
% function N = Nfun(Z,K)
%     N = ((nu*Z.*K.^alpha)./w).^(1/(1-nu));
% end       
% 
% % 3. Production function
% function Y = Yfun(Z,N,K)
%     Y = Z.*(K.^alpha).*(N.^nu);
% end  
%
% % 4. Investment
% function I = Ifun(K,Kp)
%     I = Kp - (1-delta)*K;
% end
%__________________________________________________________________________


        
end
        
        
        
