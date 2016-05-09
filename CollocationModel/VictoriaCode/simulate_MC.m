function [paths] = simulate_MC(N,Y,mt,c,v,cKS,eq,param,glob,options);

    % initial distribution over idiosync. states: draw from stat. dist
    idio_states          = gridmake(glob.pgridf,glob.agridf);
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
    
    % store policy functions
    policies    = zeros(N,options.T);
    
    % simulation will start here
    for t=1:options.T
    
        t
        % 1. define state vector
        st     = [price(:,t) prod(:,t)...
            repmat(Y(t),length(price),1) repmat(mt(t),length(price),1)];

        % 2. create basis matrix for continuation values
        glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),st(:,2));
        glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),st(:,3));
        glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),st(:,4));

        % 3. compute optimal prices and keep/change policy functions
        v  = solve_valfuncKS(c,[st(:,1) st(:,2:end)],param,glob,options);
        policies(:,t) = v.Is;
        price_opt(:,t) = v.Pc;
        
        % 4. update idiosyncratic states
        % tomorrow's prices will be the policy times inflation
        if t<options.T
            price(:,t+1)      = v.Pp.*1/exp(param.mu).*1/exp(mt(t+1));
            price_nopi(:,t+1) = v.Pp;
            log_prod     = param.rhoa.*log(prod(:,t)) + normrnd(0,param.sigmazeta,[length(price),1]);
            prod(:,t+1)  = exp(log_prod);
        end
        
    end
    
    % packup
    paths.price     = price;
    paths.prod      = prod;
    paths.policy    = policies;
    paths.price_opt = price_opt;
    paths.price_nopi= price_nopi;
        
end