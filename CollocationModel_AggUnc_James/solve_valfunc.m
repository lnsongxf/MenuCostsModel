function [v,jac] = solve_valfunc(c,s,Y,param,glob,options,xxx)
%SOLVE_VALFUNC Solves value function 
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cK,cC,cE) and the
%   stationary distribution matrix L. Solution is conditional on a
%   conjectured level of equilibrium output, Y. 
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

% cK = c(1:glob.Ns);
% cC = c(glob.Ns+1:2*glob.Ns);
% cE = c(2*glob.Ns+1:end);

cK = c(1:end/3);
cC = c(end/3+1:2*end/3);
cE = c(2*end/3+1:end);   

%__________________________________________________________________________

% Value function when keeping price
vK = valfunc('K',cE,s,[],Y,param,glob,options);

% Value function when changing price
B                       = menufun('bounds',s,[],Y,param,glob,options); 
obj                     = @(pPstar) valfunc('C',cE,s,pPstar,Y,param,glob,options);
pPstar                  = goldenx(obj,B(:,1),B(:,2));
[vC, Phi_pPAMY] = valfunc('C',cE,s,pPstar,Y,param,glob,options);


%__________________________________________________________________________
% When simulating etc, don't need the expected value function, nargin is useful
vE  = [];
if (nargin<=6)  
    % PROBLEM: How to compare vK vC to get indicator when the state space
    % is s'? Maybe can multiply by H
    
    % Expected value function
    ind = (vK > vC);    % Indicator for when value of keeping price is larger than changing price
    ind = double(ind);  % Change from logical to double
    vE = glob.Emat*(dprod(ind, glob.H*glob.Phiprime)*cK + ... 
                    dprod((1-ind), glob.H*glob.Phiprime)*cC);
end

% Compute Jacobian if requested
if (nargout==2)
    jac = [ glob.Phi,                                       zeros(glob.Ns),                                 -param.beta*glob.Phi  ;
          zeros(glob.Ns),                                   glob.Phi,                                       -param.beta*Phi_pPAMY ;
          -glob.Emat*dprod(ind, glob.H*glob.Phiprime),    -glob.Emat*dprod((1-ind), glob.H*glob.Phiprime),  glob.Phi             ];          
end

%__________________________________________________________________________
% Packup output
v.vK        = vK;
v.vC        = vC;
v.vE        = vE;
v.pPstar    = pPstar;

% v.N     = N;
% v.Y     = Y;
% v.D     = D;
% v.I     = I;
% v.AC    = AC;
% v.C     = C;




