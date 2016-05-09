function [paths] = simulate_MC(N,C,pi_sim,c,v,cKS,eq,param,glob,options);

    % initial distribution over idiosync. states: draw from stat. dist
    idio_states          = gridmake(glob.xgridf,glob.nugridf);
    s1                   = idio_states(:,1);
    s2                   = idio_states(:,2);
    state_init           = randsample(1:length(eq.L),N,true,max(eq.L,0))';
    price                = zeros(N,options.T);
    price_opt            = zeros(N,options.T);
    prod                 = zeros(N,options.T);
    price_nopi           = zeros(N,options.T);
    price(:,1)           = s1(state_init);
    price_nopi(:,1)      = s1(state_init);
    prod(:,1)            = s2(state_init);
    price_level          = zeros(1,options.T);
    price_level(1)       = exp(pi_sim(1));
    P                    = zeros(1,options.T);
    
    % store policy functions
    policies    = zeros(N,options.T);
    
    % simulation will start here
    for t=1:options.T
    
        t
        % 1. define state vector
        st     = [price(:,t) prod(:,t) repmat(C(t),length(price),1)];

        % 2. create basis matrix for continuation values
        glob.Phi_nu  = splibas(glob.nugrid0,0,glob.spliorder(2),st(:,2));
        glob.Phi_c   = splibas(glob.cgrid0,0,glob.spliorder(3),st(:,3));

        % 3. compute optimal prices and keep/change policy functions
        v  = solve_valfuncKS(c,[st(:,1) st(:,2:end)],param,glob,options);
        policies(:,t)   = v.Is;
        price_opt(:,t)  = v.Xc;
        
        % 4. update idiosyncratic states
        % tomorrow's prices will be the policy times inflation
        if t<options.T
            price(:,t+1)      = v.Xp.*1/exp(param.mu).*(1/exp(pi_sim(t+1)-param.mu)).*(1/exp(param.mu));
            price_nopi(:,t+1) = v.Xp;
            nom_prices        = (param.alpha.*param.R.*v.Xp.*price_level(t)).^(1-param.epsilon);
            P(t)              = mean(nom_prices)^(1/(1-param.epsilon));
            log_prod          = (1-param.eta).*log(prod(:,t)) + normrnd(0,param.sigmanu,[length(price),1]);
            prod(:,t+1)       = exp(log_prod);
            price_level(t+1)  = exp(log(price_level(t)) + pi_sim(t+1));
        end
        
    end
    
    % packup
    paths.price     = price;
    paths.prod      = prod;
    paths.policy    = policies;
    paths.price_opt = price_opt;
    paths.price_nopi= price_nopi;
    paths.price_level=price_level;
    paths.P         = P;
        
end