function [coeffs,R2,paths] = simulate_KS(mt,M_sim,c,v,cKS,eq,param,glob,options)

    %% initialize for simulation
    
    % distribution over idiosyncratic states
    L_sim       = zeros(length(eq.L),options.T);
    L_sim(:,1)  = eq.L;
    
    % policy functions
    pol_sim     = zeros(length(eq.L),options.T-1);
    I_sim       = zeros(length(eq.L),options.T-1);
    I_sim(:,1)  = eq.v.Is;
    
    % individual and aggregate labor and wages
    l_sim       = zeros(length(eq.L),options.T);
    ell_sim       = zeros(length(eq.L),options.T);
    Labor_sim   = zeros(1,options.T);
    w_sim       = zeros(length(eq.L), options.T);
    Wage_sim    = zeros(1,options.T);
    Ell_sim    = zeros(1,options.T);
    
    % vectors for price level, output
    P_sim       = zeros(1,options.T);
    Y_sim       = zeros(1,options.T);
    P_sim(1)    = (1/eq.Y)*M_sim(1);
    Y_sim(1)    = M_sim(1)/P_sim(1);
    
    % initialize labor and wages
    st              = gridmake(glob.pgridf,glob.agridf);
    yj              = Y_sim(1).*(eq.v.Pp).^(-param.theta);
    ellj            = (yj./st(:,2)).^(1/param.alpha);
    lj              = ellj + eq.v.Is.*param.Phi;
    l_sim(:,1)      = lj;
    Labor_sim(1)    = L_sim(:,1)'*l_sim(:,1);
    w_sim(:,1)      = (param.delta.*lj.^param.phi.*Y_sim(1).^param.sigma)./(P_sim(1));
    Wage_sim(1)     = L_sim(:,1)'*w_sim(:,1);
    
    
    
    %% simulate
    
    %t=2;
    tictic = tic;
    for t = 2:options.T 

        ylb    = 0.5*Y_sim(t-1);
        yub    = 1.5*Y_sim(t-1);
%        Yin    = (1/2)*(ylb+yub);
        Pin    = P_sim(t-1);
        Pout   = Pin;

        for itery = 1:50;
            
            Yin    = (1/2)*(ylb+yub);

            % 1. define state vector
            st     = gridmake(glob.pgridf,glob.agridf,Yin,mt(t+1));

            % 3. create basis matrix for continuation values
            glob.Phi_A  = splibas(glob.agrid0,0,glob.spliorder(2),st(:,2));
            glob.Phi_Y  = splibas(glob.Ygrid0,0,glob.spliorder(3),st(:,3));
            glob.Phi_m  = splibas(glob.mgrid0,0,glob.spliorder(4),st(:,4));

            % 4. compute real price distribution from policy functions
            v  = solve_valfuncKS(c,[st(:,1) st(:,2:end)],param,glob,options);

            % 5. Multiply by P_t guess to get nominal price distribution
            nom_prices  = v.Pp.*Pin;

            % 6. Use CES to aggregate price level
            Pout        = (L_sim(:,t-1)'*(nom_prices.^(1-param.theta)))^(1/(1-param.theta));

            % 7. Market clearing
            Yout        = M_sim(t+1)./Pout;
            
            % 8. Check convergence of Y
            down        = (Yin>Yout); 
            up          = (Yin<Yout);
            ylb         = up*Yin + down*ylb;
            yub         = up*yub + down*Yin;
            if strcmp(options.eqprint,'Y') 
                 fprintf('%2i. yin:\t%2.6f\tyout:\t%2.6f\tt:%2.1f\n',itery,Yin,Yout,toc(tictic));
            end
%            if abs(Yout-Yin)<options.toly;fprintf('Solved\n');break;end;
            if abs(Yout-Yin)<options.toly;break;end;
%             Yin = 0.9*Yin + 0.1*Yout;

        end
        
        pi              = Pout/P_sim(t-1);
        Y_sim(t)        = Yout;
        P_sim(t)        = Pout;
        pol_sim(:,t)    = v.Pp;
        I_sim(:,t)      = v.Is;
        yj              = Y_sim(t).*(v.Pp).^(-param.theta);
        ellj            = (yj./st(:,2)).^(1/param.alpha);
        lj              = ellj + v.Is.*param.Phi;
        l_sim(:,t)      = lj;
        ell_sim(:,t)    = ellj;
        
        % 9. Update distributions
        fspaceerg     = fundef({'spli',glob.pgridf,0,1});
        Pp            = max(min(v.Pp,max(glob.pgridf)),min(glob.pgridf)).*(1/pi);
        Qp            = funbas(fspaceerg,Pp);
        L_sim(:,t)    = dprod(glob.QA,Qp)'*L_sim(:,t-1);
        Labor_sim(t)  = L_sim(:,t)'*l_sim(:,t);    
        w_sim(:,t)    = (param.delta.*lj.^param.phi.*Y_sim(t).^param.sigma)./(P_sim(t));
        Wage_sim(t)   = L_sim(:,t)'*w_sim(:,t);
        Ell_sim(t)      = L_sim(:,t)'*ellj;
        
    end
    
    
    %% Finishing up...
    
    % compute regression coefficients
    X = [ones(options.T-2,1) log(Y_sim(1:options.T-2))' mt(3:options.T)'];
    y = log(Y_sim(2:options.T-1))';
    [b,bint,r,rint,stats] = regress(y,X);
    coeffs = b';
    R2     = stats(1);
    
    % packup the simulation paths
    paths.mt    = mt;
    paths.L     = L_sim;
    paths.M     = M_sim;
    paths.P     = P_sim;
    paths.Y     = Y_sim;
    paths.pol   = pol_sim;
    paths.I     = I_sim;
    paths.l     = l_sim;
    paths.Labor = Labor_sim;
    paths.w     = w_sim;
    paths.Wage  = Wage_sim;
    paths.ell   = ell_sim;
    paths.Ell   = Ell_sim;
    
end