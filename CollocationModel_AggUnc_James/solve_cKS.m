function [c,v] = solve_cKS(cold,cKS,param,glob,options)
%SOLVE_CKS Solves the value function specifically for the Krussel Smith alogrithm
%-------------------------------------------------
%
%   INPUTS
%   - cold      = Initialization of value function approximation coefficients
%   - cKS       = parameters for the law of motion (lnY = b0 + b1lnY-1 + b2 Dm)
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%
%   OUTPUTS
%   - c         = Solved value function approximation coefficients 
%   - v         = Sovled value function and its elements
%-------------------------------------------------

%% Initialize
s           = glob.s;
totaltic    = tic;

%% Compute Emat (expectations matrix on total state space)
% Simply run setup again but with different law of motion parameters.
% Do this because all of the infrastructure is inside the setup file: with
% different cKS params, we compute the new Tauchen VAR for (dM, Y), update
% the state space, construct new expectations matrices, etc. 

% fprintf('Setup\n');
% [param,glob] = setup(param,glob,cKS,options);    
% fprintf('Setup complete\n');
 

% NOTE: May need to write this out individually since the new expectations
% matrices are actually created in the calling function 'solve_KS'...

%% Compute forecast for output (In Simon's code this is prices)

DM          = s(:,3);
% X           = [ones(size(Y)),log(Y),log(DM)];
% Y           = exp(X*cKS);     % KS forecast for output

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;    
    v           = solve_valfuncKS(cold,s,param,glob,options); 
    % 2. Update c
    cK          = glob.Phi\full(v.vK);      % Note: 'full' re-fills a sparse matrix for computations
    cC          = glob.Phi\full(v.vC);
    cE          = glob.Phi\full(v.vE);    
    c           = [cK;cC;cE];
    % 3. Compute distance and update
    dc          = norm(c-cold)/norm(cold); 
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
end

%% Newton iterations

if strcmp(options.print,'Y');
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
eq.flag.cconv = false;
for citer = (1:options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfuncKS(cold,s,param,glob,options);
    % 2. Update c 
    cKold       = cold(1:glob.Ns); 
    cCold       = cold(glob.Ns+1:2*glob.Ns);
    cEold       = cold(2*glob.Ns+1:end);
    c           = cold - jac\([glob.Phi*cKold - full(v.vK) ;
                               glob.Phi*cCold - full(v.vC) ;
                               glob.Phi*cEold - full(v.vE)]);  
    % 3. Compute distances and update
    dc          = norm(c-cold)/norm(cold);
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end


