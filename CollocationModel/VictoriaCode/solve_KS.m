function [c,v] = solve_KS(cKS,param,glob,options)

    %% Unpack
    ns          = size(glob.s,1); 

    %% Initialize guesses (if val.cresult has an old guess in it, use that)
    ckold       = zeros(ns,1);  
    ccold       = zeros(ns,1); 
    ceold       = zeros(ns,1);
    cold        = [ckold;ccold;ceold];

    %% Solve value function problem
    [c,v]       = solve_cKS(cold,cKS,param,glob,options);

end