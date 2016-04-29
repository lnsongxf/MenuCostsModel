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
B                       = menufun_noagg('bounds',s,[],Y,param,glob,options); 
obj                     = @(pPstar)valfunc_noagg('C',cE,s,pPstar,Y,param,glob,options);
pPstar                  = goldenx(obj,B(:,1),B(:,2));
[vC, Phi_pPA] = valfunc_noagg('C',cE,s,pPstar,Y,param,glob,options);

pPdist = [s(:,1); pPstar];      % For Calvo, some change, some don't. Need to keep track of both sets of prices

% ind = double((vK > vC));    % Indicator for when value of keeping is larger than changing price
% pPdist = ind.*s(:,1) + (1-ind).*pPstar;    % distribution of real prices given state vector


%% Compute vE and jacobian if requested
vE  = [];
if (nargin<=6)  
    % Expected value function    
%     vE          = glob.Emat*max(vK,vC);
    s_prime     = [s(:,1)*(1/exp(param.mu)), s(:,2)];
    Phiprime    = funbas(glob.fspace,s_prime);
    vE = glob.Emat*(param.lambda*Phiprime*cK + (1-param.lambda)*Phiprime*cC);    % (dprod(ind, Phiprime)*cK + dprod((1-ind), Phiprime)*cC);
end

if (nargout==2)
    jac = [ glob.Phi,                               zeros(glob.Ns),          -param.beta*glob.Phi;
          zeros(glob.Ns),                           glob.Phi,                   -param.beta*Phi_pPA ;
          -glob.Emat*param.lambda*glob.Phiprime,    -glob.Emat*(1-param.lambda)*glob.Phiprime,   glob.Phi       ];          
end

%% Packup output

% Value function coefficients
v.vK        = vK;
v.vC        = vC;
v.vE        = vE;

% Optimal policy functions (correspond to each state element in 's').
v.pPstar    = pPstar;   % optimal price if changing at given state
v.pPdist    = pPdist;   % distribution of prices across changers and non-changers
% v.ind       = ind;      % Who did/didn't change prices
v.ystar     = menufun_noagg('output',s,pPdist,Y,param,glob,options);  % Use dist, because want values at *actual* prices
v.nstar     = menufun_noagg('labour',s,pPdist,Y,param,glob,options);
% v.wPstar    = menufun_noagg('realwage',s,pPdist,Y,param,glob,options);





