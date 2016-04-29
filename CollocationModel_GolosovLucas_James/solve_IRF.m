function [ xxx ] = solve_IRF(eq,param,glob,options)
%SOLVE_STATDIST Compute the stationary distribution and implied output
%----------------------------------------------------------------
%   INPUTS
%
%
%   OUTPUTS
%----------------------------------------------------------------

%
% NOTE: Big problem: for some reasson feeding the solve_noagg functijon the
% old coefficient vector does NOT return the same coefficieints. This
% causes the "equilibrium" output to be much lower than in the stationary
% state. 
%


%%
% Guess initial sequence for equilibrium output: 
% Y = Ystationary at begginning and end of path, guess basic step function
% for the rest of the path
Yguess           = zeros(options.numtrans,1);
Yguess(1)        = eq.Y;
Yguess(2:options.numtrans/2) = eq.Y; % + 0.01;
Yguess(options.numtrans/2+1:end-1) = eq.Y; % + 0.005;
Yguess(end)      = eq.Y;

% Pre-fill actual output path
Yactual         = zeros(options.numtrans,1);
Yactual(1)      = eq.Y;
Yactual(end)    = eq.Y;

% Initialize the i+1 cE coefficient vector
coeffs          = zeros(length(eq.c),options.numtrans+1);
coeffs(:,end)   = eq.c;
glob.CEi1       = coeffs(2*end/3+1:end,end);      % Set final expectations vector

% Initialize the distributions
dist            = zeros(length(eq.L),options.numtrans);
dist(:,1)       = eq.L;


% Path for wage inflation, one stddev shock in period 2, with possible
% AR(1) runoff
rhom            = 0;
piw             = param.mu*ones(options.numtrans+1,1);
shock           = zeros(options.numtrans+1,1);
shock(2)        = 0;   % 2*param.sigm;       % one stddev shock
for t = 2:options.numtrans
    piw(t)      = (1-rhom)*param.mu + shock(t);    % + rhom*piw(t-1)  
end
piw             = piw(1:end-1);



%% Solve problem recursively given piw and Ypath 
dY              = 1;
options.IRF    = 'Y';
while dY > 1e-6

    %% First iterate recursively/backwards
    for t = 1:options.numtrans-1
        glob.piw    = piw(options.numtrans - t);
        
        v           = ...
            solve_valfunc_noagg(eq.c,glob.sf,Yguess(options.numtrans-t),param,glob,options);
        % 2. Update c
        
        
        cK          = glob.Phif\full(v.vK);      % Note: 'full' re-fills a sparse matrix for computations
        cC          = glob.Phif\full(v.vC);
        
        s_primef     = [glob.sf(:,1)*(1/exp(glob.piw)), glob.sf(:,2)];
        Phiprime    = funbas(glob.fspace,s_primef);
        glob.Emat   = kron(glob.P,speye(glob.nf(1)));
        vE          = glob.Emat*(dprod(v.ind, Phiprime)*cK + dprod((1-v.ind), Phiprime)*cC);
%         vE          = glob.Emat*max(v.vK,v.vC);
        cE          = glob.Phif\full(vE);
        
        coeffs(:,options.numtrans-t)           = [cK;cC;cE]; 
        % Get t+1 expectations coefficient vector
        glob.CEi1   = coeffs(2*end/3+1:end,options.numtrans-t);  %  eq1.c(2*end/3+1:end); % eq.c(2*end/3+1:end);        
        
%         eq1          = solve_cL(Yguess(options.numtrans - t),param,glob,options);        

%         coeffs(:,options.numtrans-t) = eq1.c;
%         glob.CEi1   = coeffs(2*end/3+1:end,options.numtrans-t);  %  eq1.c(2*end/3+1:end); % eq.c(2*end/3+1:end);        
    end
    
    
    %% Second, solve forwards
    for t = 1:options.numtrans-1
        glob.piw    = piw(t);
        c           = coeffs(:,t);                % Valfun coeffs at t
        glob.CEi1   = coeffs(2*end/3+1:end,t+1); % Valfun expectation coeffs at t+1
        Yactual(t)  = Yguess(t);                 % Initialize Yactual at guess, Ypath    
        
        % Solve partial equilibrium
        for iterY = 1:50;
            Yin     = Yactual(t);
            v       = solve_valfunc_noagg(c,glob.sf,Yin,param,glob,options,1);
            Ynew = ...
                (param.alpha^(1-param.epsilon)*(v.pPdist'.^(1-param.epsilon)*dist(:,t)))^...
                (1/(param.gamma*(param.epsilon-1)));
            Yout        = options.damp*Yin   + (1-options.damp)*Ynew;
            Yactual(t)	= Yout;
            if strcmp(options.eqprint,'Y')
                fprintf('%2i. Yin:\t%2.4f\t Yout:\t%2.4f\n',iterY,Yin,Yout);
            end
            if abs(Yin-Yout)<options.tolY
                fprintf('Solved for output, t = %1i\n', t)
                break
            end
        end

        
        
        % 4. Compute next period distribution over states
        pPdist      = max(min(v.pPdist,max(glob.pPgridf)), min(glob.pPgridf)).*(1/exp(glob.piw));
        fspaceergpP = fundef({'spli',glob.pPgridf,0,1});   % Linear interpolant
        QpP         = funbas(fspaceergpP,pPdist);
        dist(:,t+1) = dprod(glob.QV,QpP)'*dist(:,t);
    end
    
    % Plot adjustment of Yactual to the guess Ypath
    figure(113)
    plot((Yguess(1:end)));    %-eq.Y)/eq.Y*100)
    hold all
    plot((Yactual(1:end)));     % -eq.Y)/eq.Y*100) 
    
%     line([0 options.numtrans],[0 0],'color','k')
    hold off 
    legend('Y guess','Y actual')
    ylabel('Percent deviation from stationary state output')
%     ylim([eq.Y-0.001 eq.Y + 0.001])
    dY              = norm(Yactual-Yguess)/norm(Yguess);    
    drawnow
    fprintf('Difference this iteration  = %f \n', dY)
  
    Yguess = options.damp*Yguess + (1 - options.damp)*Yactual;    
%     Yguess(2:end-1) = Yactual(2:end-1);     % options.damp*Ypath + (1 - options.damp)*
    
end


xxx = 1;


%% 

end





