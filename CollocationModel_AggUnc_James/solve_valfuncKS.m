function [v,jac] = solve_valfuncKS(c,s,param,glob,options,xxx)
%SOLVE_VALFUNCKS Solves value function inside the KS step
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cK,cC,cE) inside the 
%   Krussel Smith step. Solution is conditional on an implied level of 
%   next period's equilibrium output given by the current guess for the law
%   of motion parametetrs: (lnY = b0 + b1lnY-1 + b2 Dm)
%
%   INPUTS
%   - c         = current collocation coefficient matrix
%   - s         = state space
%   - Y         = conjectured value of output, Y
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = 
%   - jac       = Jacobian of the value functions  
%-------------------------------------------------

%% Unpack coefficient vector
cK = c(1:end/3);
cC = c(end/3+1:2*end/3);
cE = c(2*end/3+1:end); 

%__________________________________________________________________________
%% Solve problem 

% Value function when keeping price
vK = valfuncKS('K',cE,s,[],param,glob,options);

% Value function when changing price
B                       = menufun('bounds',s,[],[],param,glob,options); 
obj                     = @(pPstar) valfuncKS('C',cE,s,pPstar,param,glob,options);
pPstar                  = goldenx(obj,B(:,1),B(:,2));
[vC, Phi_pPAMY] = valfuncKS('C',cE,s,pPstar,param,glob,options);

ind = (vK > vC);    % Indicator for when value of keeping price is larger than changing price
ind = double(ind);  % Change from logical to double
v.ind       = ind;  % Record who did/didn't change prices

pPdist = ind.*s(:,1) + (1-ind).*pPstar;    % distribution of real prices given state vector

% Compute vE and jacobian if requested
vE  = [];
if (nargin<=5)  
    % Expected value function
    
    ind = kron(ind, ones(glob.Ny*glob.Ny*glob.Nm,1));       % Match dimensions of Phiprime
    
    vE = glob.Emat*glob.H*(dprod(ind, glob.Phiprime)*cK + ... 
                    dprod((1-ind), glob.Phiprime)*cC);
end
if (nargout==2)
    jac = [ glob.Phi,                                       zeros(glob.Ns),                                 -param.beta*glob.Phi  ;
          zeros(glob.Ns),                                   glob.Phi,                                       -param.beta*Phi_pPAMY ;
          -glob.Emat*glob.H*dprod(ind, glob.Phiprime),    -glob.Emat*glob.H*dprod((1-ind), glob.Phiprime),  glob.Phi             ];          
end

%% Packup output

% Value function coefficients
v.vK        = vK;
v.vC        = vC;
v.vE        = vE;

% Optimal policy functions (correspond to each state element in 's').
v.pPstar    = pPstar;   % optimal price if changing at given state
v.pPdist    = pPdist;   % distribution of prices across changers and non-changers
v.ystar     = menufun('output',s,pPdist,v.ind,param,glob,options);  % Use dist, because want values at *actual* prices, not optimal if they were changing
v.lstar     = menufun('labour',s,pPdist,v.ind,param,glob,options);
v.wPstar    = menufun('realwage',s,pPdist,v.ind,param,glob,options);




