function [c,v] = solve_cKS(cold,cKS,param,glob,options)

s           = glob.s;
totaltic    = tic;

%% Compute Emat (expectations matrix)
K           = glob.s(:,2);
A           = glob.s(:,4);
X           = [ones(size(K)),log(K),log(A),log(A).*log(K)];
Kp          = exp(X*cKS.gamma_K);
Kp          = max(min(Kp,max(glob.Kgrid)),min(glob.Kgrid));
Phi_K       = splibas(glob.Kgrid0,0,glob.spliorder(2),kron(Kp,glob.iNe));  
Phi_KzA     = dprod(glob.Emat_Phi_zA,Phi_K);
Phi_kKzA    = dprod(Phi_KzA,glob.Emat_Phi_k);
Ikronw      = glob.Emat_Ikronw;
size(Phi_K)
glob.Emat   = Ikronw*Phi_kKzA;  

%% Compute prices
K           = s(:,2); 
A           = s(:,4);
X           = [ones(size(K)),log(K),log(A),log(A).*log(K)];
P           = exp(X*cKS.gamma_p);     % KS forecast
% P           = glob.p*ones(size(K));   %%%% DEBUG

%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfuncKS(cold,s,P,param,glob,options); 
    % 2. Update c
    c1          = glob.Phi\full(v.v1);   
    c2          = glob.Phi\full(v.v2); 
    c           = [c1;c2];
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
    [v,jac]     = solve_valfuncKS(cold,s,P,param,glob,options);  
    % 2. Update c 
    c1old       = cold(1:glob.Ns); 
    c2old       = cold(glob.Ns+1:end);
    c           = cold - jac\([glob.Phi*c1old - full(v.v1) ; 
                               glob.Phi*c2old - full(v.v2)]);   
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
