function eq  = solve_xL(cbar,param,glob,options)

%% A. Globals 
s           = glob.s;  
sf          = glob.sf;
xgrid       = glob.xgrid;
ns          = size(s,1);

%% B. Compute equilibrium objects that depend on cbar
% ----- None -----

%% Initialise guesses (if val.cresult has an old guess in it, use that)
ckold       = zeros(ns,1);
ccold       = zeros(ns,1);
ceold       = zeros(ns,1);
c           = options.cresult;
if ~isempty(c) && strcmp(options.Loadc,'Y')
    ckold    = c(1:end/3);
    ccold    = c(end/3+1:2*end/3);
    ceold    = c(2*end/3+1:end);
end
cold        = [ckold;ccold;ceold];
totaltic    = tic;

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc_GL(cold,s,cbar,param,glob,options); 
    % 2. Update c
    ck          = funbas(glob.fspace,s)\full(v.vk);   
    cc          = funbas(glob.fspace,s)\full(v.vc); 
    ce          = funbas(glob.fspace,s)\full(v.ve); 
    c           = [ck;cc;ce]; 
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
    [v,jac]     = solve_valfunc_GL(cold,s,cbar,param,glob,options);
    % 2. Update c 
    ckold = cold(1:end/3);  
    ccold = cold(end/3+1:2*end/3);
    ceold = cold(2*end/3+1:end);
    c     = cold - jac\([funbas(glob.fspace,s)*ckold - full(v.vk) ; 
          	funbas(glob.fspace,s)*ccold - full(v.vc) ;
         	funbas(glob.fspace,s)*ceold - full(v.ve)]);   
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
    if eq.flag.cconv,break,end;
end
    
%% Solve again on a finer grid for p
%  NOTE (VG): for the stationary distribution
glob.Phi_nu      = glob.Phi_nuf; 
v                = solve_valfunc_GL(c,sf,cbar,param,glob,options);


%% Compute stationary distribution
%  NOTE (VG): I think he's avoiding policies off the max of the grid
Xp              = min(v.Xp,max(xgrid));     
fspaceergx      = fundef({'spli',glob.xgridf,0,1});
QX              = funbas(fspaceergx,Xp); 
QNu             = glob.QNu;
Q               = dprod(QNu,QX); 

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

% Pa = (L'*((v.Pp).^(1-param.theta))).^(1/(1-param.theta));
% Ca = 1/Pa;

cbara = (param.alpha.^(1-param.epsilon).*(sf(:,1)'.^(1-param.epsilon)*L))...
    .^(1/(param.gamma*(param.epsilon - 1)));

%% Pack-up output
eq.v        = v;
eq.c        = c;
eq.L        = L;
eq.Q        = Q;
eq.cbar     = cbara;

end