function eq  = solve_pL(Y,param,glob,options)

%% A. Globals 
s           = glob.s;  
sf          = glob.sf;
pgrid       = glob.pgrid;
ns          = size(s,1);

%% B. Compute equilibrium objects that depend on Y
% ----- None -----

%% Initialize guesses (if val.cresult has an old guess in it, use that)
%  Here, we're going to solve for the expected value function only
ceold       = zeros(ns,1);
c           = options.cresult;
if ~isempty(c) && strcmp(options.Loadc,'Y')
    ceold    = c(1:ns);
end
totaltic    = tic;

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc_menucost(ceold,s,Y,param,glob,options); 
    % 2. Update c
    c          = glob.Phi\full(v.vf); 
    % 3. Compute distance and update
    dc          = norm(c-ceold)/norm(ceold); 
    ceold        = c;
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
    [v,jac]     = solve_valfunc_menucost(ceold,s,Y,param,glob,options);
    % 2. Update c 
    c           = ceold - jac\(glob.Phi*ceold - full(v.vf));  
    % 3. Compute distances and update
    dc          = norm(c-ceold)/norm(ceold);
    ceold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv,break,end;
end
    
%% Solve again on a finer grid for p
%  NOTE (VG): for the stationary distribution
glob.Phi_A      = glob.Phi_Af; 
v               = solve_valfunc_menucost(c,sf,Y,param,glob,options);


%% Compute stationary distribution
%  NOTE (VG): I think he's avoiding policies off the max of the grid
Pp              = min(v.Pp,max(pgrid));     
fspaceergp      = fundef({'spli',glob.pgridf,0,1});
QP              = funbas(fspaceergp,Pp); 
QA              = glob.QA;
Q               = dprod(QA,QP); 

% [vv,dd]         = eigs(Q');
% dd              = diag(dd);
% Lv              = vv(:,dd==max(dd));
% L               = Lv/sum(Lv);
L               = ones(size(Q,1),1);
L               = L/sum(L);

for itL = (1:options.itermaxL);
    Lnew    = Q'*L;  
    dL      = norm(Lnew-L)/norm(L);  
    if (dL<options.tolL),break,end;
    if mod(itL,100)==0 
        if strcmp(options.print,'Y')
            fprintf('dL:\t%1.3e\n',dL);
        end
    end
    L       = Lnew;
end

%% Compute aggregates

Pa = (L'*((v.Pp).^(1-param.theta))).^(1/(1-param.theta));
Ca = 1/Pa;

%% Pack-up output
eq.v    = v;
eq.c    = c;
eq.L    = L;
eq.Pa   = Pa;
eq.Ca   = Ca;
eq.Q    = Q;
eq.Y    = Ca;

end