function out = menufun(flag,s,pPstar,ind,Y,param,glob,options)
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
gamma   = param.gamma;          % risk aversion
epsilon = param.epsilon;        % elasticity of substition
alpha   = param.alpha;          % disutility of labour
rhov    = param.rhov;           % AR(1) coefficient for productivity
sigv    = param.sigv;           % Std dev of productivity process
k       = param.k;              % menu cost 
mu      = param.mu;             % Drift parameter for money growth
sigm    = param.sigm;           % Std dev of money growth shocks
R       = param.R;              % Stationary interest rate


%% Cases
switch flag
    case 'bounds'
        pPpmin   = ones(size(s,1),1)*min(glob.pPgrid); 
        pPpmax   = ones(size(s,1),1)*max(glob.pPgrid);
        out     = [pPpmin,pPpmax];
    case 'PIK'
        pP      = s(:,1);
        V       = s(:,2);
        out = Y.^(1 - epsilon*gamma).*(alpha*pP).^(-epsilon).*(pP - 1./V);
    case 'PIC'
        V       = s(:,2);
        MenuCost = MCfun();
        out = Y.^(1 - epsilon*gamma).*(alpha*pPstar).^(-epsilon).*(pPstar - 1./V) - MenuCost;
%     case 'output'
%         out       = Yfun(pPstar,Y);         % Firm's output 
%     case 'labour'
%         V         = s(:,2);
%         ystar     = Yfun(pPstar,Y);
%         out       = Nfun(ystar, V, ind);    % Firm's labour demand
%     case 'realwage'
%         V         = s(:,2);
%         ystar     = Yfun(pPstar,Y);
%         nstar     = Nfun(ystar, V, ind);       
%         out       = Wfun(nstar,Y);          % Firm's real wage        
end


%% NESTED FUNCTIONS

% 1. Menu costs
function MenuCost = MCfun()
    switch options.MC
        case 'Y'
            MenuCost = k*ones(size(s,1),1); 
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
        
        
        
