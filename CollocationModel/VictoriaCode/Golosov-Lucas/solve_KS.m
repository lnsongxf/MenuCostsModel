function [c,v,KS_coeffs,R2,paths] = solve_KS(cKS,eq,param,glob,options)

    %% Setup problem (this stuff has to change when coefficients change)
    fprintf('Setup\n');
    [param,glob]    = setup_ks(cKS,param,glob,options);      
    fprintf('Setup complete\n');

    %% Unpack
    ns          = size(glob.s,1); 

    %% Initialize guesses (if val.cresult has an old guess in it, use that)
    ckold       = zeros(ns,1);  
    ccold       = zeros(ns,1); 
    ceold       = zeros(ns,1);
    cold        = [ckold;ccold;ceold];

    %% Solve value function problem
    [c,v]       = solve_cKS(cold,cKS,param,glob,options);
    
    %% Set-up for simulations
    
    % first draw inflation shocks: discrete
    rng(219);
    %pi_sim      = param.mu + randsample(glob.supp_e,options.T,true,glob.f);
    pi_sim       = param.mu + normrnd(0,param.sigmam,[options.T,1]);
    
    % simulate and get coefficients
    [KS_coeffs, R2, paths] = simulate_KS(pi_sim,c,v,cKS,eq,param,glob,options);

    
end