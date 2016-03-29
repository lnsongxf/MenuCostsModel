function out = menufun_noagg(flag,s,pPstar,ind,Y,param,glob,options)
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
    case 'output'
        out       = Yfun(pPstar,Y);         % Firm's output 
    case 'labour'
        A         = s(:,2);
        ystar     = Yfun(pPstar,Y);
        out       = Nfun(ystar, A, ind);    % Firm's labour demand
    case 'realwage'
        A         = s(:,2);
        ystar     = Yfun(pPstar,Y);
        nstar     = Nfun(ystar, A, ind);       
        out       = Wfun(nstar,Y);          % Firm's real wage        
end


%% NESTED FUNCTIONS

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
function nstar = Nfun(ystar,A,ind)
    nstar = (ystar./A).^(1/alpha) + ind*Phicost; 
end       

% 3. Production function
function ystar = Yfun(pPstar,Y)
    ystar = Y.*(pPstar).^(-theta);
end

% 4. Real wage paid by firm
function wPstar = Wfun(lstar,Y)
    wPstar = delta*lstar.^(phielas).*Y.^(sigma);
end


        
end
        
        
        
