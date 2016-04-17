function [coeffs,R2,paths] = simulate_KS(pi_sim,c,v,cKS,eq,param,glob,options)

    %% initialization
    
    % distribution over idiosyncratic states
    L_sim       = zeros(length(eq.L),options.T);
    L_sim(:,1)  = eq.L;

    % policy functions
    pol_sim     = zeros(length(eq.L),options.T-1);
    I_sim       = zeros(length(eq.L),options.T-1);
    
    % aggregate variable: C
    C_sim       = zeros(1,options.T);
    
    %% simulation
    
    %t = 1;
    tictic = tic;
    for t = 1:options.T
        
        % bisect over aggregate variable C
        clb    = 0.5*eq.cbar;
        cub    = 1.5*eq.cbar;
        
        for iterc = 1:50;
            
            cin    = (1/2)*(clb+cub);
            
            % 1. define state vector
            st     = gridmake(glob.xgridf,glob.nugridf,cin);
            
            % 2. create basis matrices for continuation values
            glob.Phi_nu  = splibas(glob.nugrid0,0,glob.spliorder(2),st(:,2));
            glob.Phi_c   = splibas(glob.cgrid0,0,glob.spliorder(3),st(:,3));
            
            % 3. compute x distribution from policy functions
            v  = solve_valfuncKS(c,[st(:,1) st(:,2:end)],param,glob,options);
            
            % 4. compute implied C
            supp = (cin.^(1-param.epsilon*param.gamma).*...
              (param.alpha.*v.Xp).^(-param.epsilon)).^(1-1./param.epsilon);
            cout = (supp'*L_sim(:,t)).^(param.epsilon./(param.epsilon - 1));
            
            % 5. Check convergence of C
            down        = (cin>cout); 
            up          = (cin<cout);
            clb         = up*cin + down*clb;
            cub         = up*cub + down*cin;
            if strcmp(options.eqprint,'Y') 
                 fprintf('%2i. cin:\t%2.6f\tcout:\t%2.6f\tt:%2.1f\n',iterc,cin,cout,toc(tictic));
            end
            if abs(cout-cin)<options.tolcagg;break;end;

        end
        
        % update C
        C_sim(t) = cout;
        
        % update distribution over individual states
        % construct Q_x
%         if t < options.T
%             fspaceerg     = fundef({'spli',glob.xgridf,0,1});
%             Q_X           = sparse(length(eq.L),length(glob.xgridf));
%             for e=1:glob.Ne
%                 Xp        = max(min(v.Xp,max(glob.xgridf)),min(glob.xgridf)).*(1/exp(glob.supp_e(e))).*(1/exp(param.mu));
%                 Q_Xe      = funbas(fspaceerg,Xp); 
%                 Q_X       = Q_X + glob.f(e).*Q_Xe;
%             end
%             L_sim(:,t+1)  = dprod(glob.QNu,Q_X)'*L_sim(:,t);
%         end
        
        if t < options.T
            
            fspaceerg     = fundef({'spli',glob.xgridf,0,1});
            
            % find what shock we're going to tomorrow
            ind       = find(pi_sim(t+1) == glob.supp_e+param.mu);
            Xp        = max(min(v.Xp,max(glob.xgridf)),min(glob.xgridf)).*...
                        (1/exp(glob.supp_e(ind))).*(1/exp(param.mu));
            Q_X           = funbas(fspaceerg,Xp); 
            L_sim(:,t+1)  = dprod(glob.QNu,Q_X)'*L_sim(:,t);
            
        end
        
    end

    
    %% Finishing up...
    
    % compute regression coefficients
    X = [ones(options.T-1 - options.Tburn+1,1) log(C_sim(options.Tburn:options.T-1))' pi_sim(options.Tburn+1:end)];
    y = log(C_sim(options.Tburn+1:end))';
    [b ,bint,r,rint,stats] = regress(y,X);
    coeffs = b';
    R2     = stats(1);
    
    % packup the simulation paths
    paths.logpi    = pi_sim(2:end);
    paths.L        = L_sim;
    paths.C        = C_sim;

end