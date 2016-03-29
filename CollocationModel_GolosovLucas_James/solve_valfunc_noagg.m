function [v,jac] = solve_valfunc_noagg(c,s,Y,param,glob,options,xxx)
%SOLVE_VALFUNC_NOAGG Solves value function in no uncertainty case
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

%% Unpack coefficient vector
cK = c(1:end/3);
cC = c(end/3+1:2*end/3);
cE = c(2*end/3+1:end);   

%% Value function when keeping price
vK = valfunc_noagg('K',cE,s,[],Y,param,glob,options);

%% Value function when changing price
B                       = menufun('bounds',s,[],[],Y,param,glob,options); 
obj                     = @(pPstar)valfunc_noagg('C',cE,s,pPstar,Y,param,glob,options);
pPstar                  = goldenx(obj,B(:,1),B(:,2));
[vC, Phi_pPV] = valfunc_noagg('C',cE,s,pPstar,Y,param,glob,options);

ind = (vK > vC);    % Indicator for when value of keeping price is larger than changing price
ind = double(ind);

pP = s(:,1);
pPdist = ind.*pP + (1-ind).*pPstar;    % distribution of real prices given state vector


%% Compute vE and jacobian if requested
vE  = [];
if (nargin<=6)  
    % Expected value function    
    vE = glob.Emat*(dprod(ind, glob.Phiprime)*cK + dprod((1-ind), glob.Phiprime)*cC);
end

if (nargout==2)
    jac = [ glob.Phi,                               zeros(glob.Ns),          -param.beta*glob.Phi;
          zeros(glob.Ns),                        glob.Phi,                   -param.beta*Phi_pPV ;
          -glob.Emat*dprod(ind, glob.Phiprime),    -glob.Emat*dprod((1-ind), glob.Phiprime),  glob.Phi       ];          
end

%% Packup output

% Value function coefficients
v.vK        = vK;
v.vC        = vC;
v.vE        = vE;

% Optimal policy functions (correspond to each state element in 's').
v.pPstar    = pPstar;   % optimal price if changing at given state
v.pPdist    = pPdist;   % distribution of prices across changers and non-changers
v.ind       = ind;      % Who did/didn't change prices
% v.ystar     = menufun('output',s,pPdist,ind,Y,param,glob,options);  % Use dist, because want values at *actual* prices, not optimal if they were changing
% v.nstar     = menufun('labour',s,pPdist,ind,Y,param,glob,options);
% v.wPstar    = menufun('realwage',s,pPdist,ind,Y,param,glob,options);





