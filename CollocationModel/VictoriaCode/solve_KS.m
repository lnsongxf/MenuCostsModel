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
    
    % where the coefficients will be saved from each simulation
    KS_coeffs      = zeros(options.S,numel(fieldnames(cKS)));
    rng(219);
    for s=1:options.S
        
        % draw money shocks
        mt          = zeros(1,options.T+1);
        mt(1)       = param.mu;
        for t = 2:options.T+1;
            mt(t) = param.mu*(1-param.rhom) + param.rhom*mt(t-1) + param.sigmaeps*randn;
            mt(t) = max(min(mt(t),max(glob.mgrid)),min(glob.mgrid));
        end 
        Minit       = 0;
        M_sim       = zeros(1,options.T+1);
        M_sim(1)    = Minit + mt(1);
        for t = 2:options.T+1
            M_sim(t) = mt(t) + M_sim(t-1);
        end
        M_sim       = exp(M_sim);
        
        % simulate and compute regression coefficients
        [KS_coeffs(s,:),R2,paths] = simulate_KS(mt,M_sim,c,v,cKS,eq,param,glob,options);
        s
    end
    
    

end