function [param,glob] = setup_ks(cKS,param,glob,options)

%% State space for a, Y, m: set up all at once

Np = glob.n(1);
Na = glob.n(2);
Ny = glob.n(3);
Nm = glob.n(4);

%------------------------------

% James' way: COMPARE TO MY WAY
[Pa,agrid,Pssa]      = setup_MarkovZ(Na,param.sigmazeta,param.rhoa,1);
muvar = [cKS.b0 + cKS.b2*param.mu*(1-param.rhom); param.mu*(1-param.rhom)];
Avar = [param.rhom, 0; cKS.b2*param.rhom, cKS.b1]';
Svar = [cKS.b2*param.sigmaeps; param.sigmaeps];
[tmpgrid,Pym,Pssym]              = tauchenvar([Ny; Nm],muvar,Avar,Svar);
Pym = Pym';     % Make sure 
mgrid = exp(unique(tmpgrid(2,:)))'; % Take D(ln(M_t)) back to DM_t 
Ygrid = exp(unique(tmpgrid(1,:)))'; % Take ln(Y_t) back to Y_t    
agrid0      = agrid;
Ygrid0      = Ygrid;
mgrid0      = mgrid;

%------------------------------

% % My way
% % Folding a, Y, m into a single VAR
% Nvar = [Na; Ny; Nm];
% muvar = [0; ...
%     cKS.b0 + cKS.b2*param.mu*(1-param.rhom); ...
%     param.mu*(1-param.rhom)];
% Avar = [param.rhoa, 0, 0; ...
%     0, cKS.b1, cKS.b2*param.rhom; ...
%     0, 0, param.rhom];
% Svar = [param.sigmazeta; ...
%     cKS.b2*param.sigmaeps; ...
%     param.sigmaeps];
% 
% % Create grid and transition matrix for a, Y, m
% [aYm_grid,A,Assa] = tauchenvar(Nvar,muvar,Avar,Svar);
% aYm_grid          = exp(aYm_grid');
% A                 = A';
% Assa              = Assa';
% emat_v          = kron(A,speye(Np));

% Isolate separate a, Y, m grids
% agrid       = unique(aYm_grid(:,1));
% agrid0      = agrid;
% Ygrid       = unique(aYm_grid(:,2));
% Ygrid0      = Ygrid;
% mgrid       = unique(aYm_grid(:,3));
% mgrid0      = mgrid;

%% State space for endogenous variable p
Np              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
pgrid           = nodeunif(Np,glob.pmin.^curv(1),glob.pmax.^curv(1)).^(1/curv(1));  % Adds curvature
pgrid0          = pgrid;    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',pgrid,0,glob.spliorder(1)},...
                         {'spli',agrid,0,glob.spliorder(2)},...
                         {'spli',Ygrid,0,glob.spliorder(3)},...
                         {'spli',mgrid,0,glob.spliorder(4)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
%pgrid           = s(s(:,2)==s(1,2),1); 
%agrid           = s(s(:,1)==s(1,1),2);
pgrid           = unique(s(:,1));
agrid           = unique(s(:,2));
Ygrid           = unique(s(:,3));
mgrid           = unique(s(:,4));
Np              = size(pgrid,1); 
Na              = size(agrid,1);
Ny              = size(Ygrid,1); 
Nm              = size(mgrid,1);

%% Create inflation grid and expectations matrix (E in my notes)

% inflation grid
yypmp       = gridmake(Ygrid, Ygrid, mgrid);
pigrid      = yypmp(:,3) - log(yypmp(:,2)) + log(yypmp(:,1));
Npi         = Ny*Ny*Nm;

% Form the matrix that determines whether you can transit from:
% Y --> Y', m', pi'
Pryypmppp   = sparse(Ny,Ny*Nm*Npi);
p = 1;
for m = 1:Nm*Ny
    for y = 1:Ny
        temp = sparse(Ny,Ny*Nm);
        temp(y,m) = 1;
        Pryypmppp(:,p:p+Ny*Nm-1) = temp;
        p=p+Ny*Nm;
    end
end

% First Kroenecker product: repeat vertically across m dimension
Prymypmppp     = kron(Pryypmppp,ones(Nm,1));
% Second Kroenecker product: repeat horizontally and vertically for a dim
Praymypmppp    = kron(ones(Na),Prymypmppp);
% Third Kroenecker product: repeat horizontally and vertically for p dim
Prpaympypmppp  = kron(ones(Np),Praymypmppp);

% Normalize to reflect that each transition is equally probable
Prpaympypmppp = sparse(diag(1./sum(Prpaympypmppp,2))*Prpaympypmppp);

% This is the matrix for transitioning between the N states
%N_trans = kron(A,speye(Np));
N_trans = kron(kron(Pym,Pa),speye(Np));

% I *think* this is the matrix that should be used to form expectations
glob.exp_matrix = sparse(N_trans*Prpaympypmppp);

% %% Compute expectations matrix
% 
% Phi             = funbas(fspace,s);
% glob.Emat       = kron(A,speye(Np))*Phi;

%% Construct fine grid for histogram

% % price grid
% pgridf          = nodeunif(glob.nf(1),glob.pmin.^glob.curv(1),glob.pmax.^glob.curv(1)).^(1/glob.curv(1));
% Npf             = size(pgridf,1);
% 
% % exogenous states
% Naf                 = glob.nf(2);
% Nyf                 = glob.nf(3);
% Nmf                 = glob.nf(4);
% Nvarf               = [Naf; Nyf; Nmf];
% [aYm_gridf,~,~]     = tauchenvar(Nvarf,muvar,Avar,Svar);
% aYm_gridf           = exp(aYm_gridf');
% agridf              = unique(aYm_gridf(:,1));
% Ygridf              = unique(aYm_gridf(:,2));
% mgridf              = unique(aYm_gridf(:,3));
% % agridf          = agrid;     
% % Naf             = size(agridf,1);
% 
% % whole state space
% sf              = gridmake(pgridf,agridf,Ygridf,mgridf);
% Nsf             = size(sf,1);
% 
% % store the fine grids
% glob.pgridf     = pgridf;
% glob.agridf     = agridf;
% glob.Ygridf     = Ygridf;
% glob.mgridf     = mgridf;
% glob.sf         = sf;
% glob.Nsf        = Nsf;

%% Compute QA matrix for approximation of stationary distribution
% glob.QA         = kron(A,ones(Npf,1)); 

%% Create one time only basis matrices: these might need to be changed
glob.Phi_A      = splibas(agrid0,0,spliorder(2),s(:,2));                % Used in Bellman / Newton computing expected values
% glob.Phi_Af     = splibas(agrid0,0,spliorder(2),sf(:,2));               % Used when solving on fine grid
glob.Phi_Y      = splibas(Ygrid0,0,spliorder(2),s(:,3));                % Used in Bellman / Newton computing expected values
% glob.Phi_Yf     = splibas(Ygrid0,0,spliorder(2),sf(:,3));               % Used when solving on fine grid
glob.Phi_m      = splibas(mgrid0,0,spliorder(2),s(:,4));                % Used in Bellman / Newton computing expected values
% glob.Phi_mf     = splibas(mgrid0,0,spliorder(2),sf(:,4));               % Used when solving on fine grid
Phi_P           = splibas(pgrid0,0,spliorder(1),s(:,1));

% Phi(s):
glob.Phi        = dprod(glob.Phi_m,dprod(glob.Phi_Y,dprod(glob.Phi_A,Phi_P))); 
% Phi(stilde):
pterm           = kron(s(:,1),ones(Npi,1));
piterm          = kron(ones(Np*Na*Ny*Nm,1),1./pigrid);
ppiterm         = pterm.*piterm;
sprime          = kron(s(:,2:end),ones(Npi,1));
%ppiterm         = kron(kron(kron(s(:,1), Ygrid.^(-1)), Ygrid), mgrid.^(-1));
glob.Phi_stilde = funbas(fspace,[ppiterm sprime]);

%% Declare additional global variables
glob.pgrid0     = pgrid0;
glob.pgrid      = pgrid;
glob.agrid0     = agrid0;
glob.agrid      = agrid;
glob.Ygrid0     = Ygrid0;
glob.Ygrid      = Ygrid;
glob.mgrid0     = mgrid0;
glob.mgrid      = mgrid;
% glob.A          = A;
% glob.Assa       = Assa;
glob.Np         = Np;
glob.Na         = Na;
glob.Ny         = Ny;
glob.Nm         = Nm;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;
glob.Npi        = Npi;

end