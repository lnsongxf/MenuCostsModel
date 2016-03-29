function [param,glob] = setup_ss(param,glob,options)

%% State space for idiosyncratic productivity
% One persistent shock
Na              = glob.n(2); 
[A,agrid,Assa]  = setup_MarkovZ(Na,param.sigmazeta,param.rhoa,1);         
agrid0          = agrid;

%% State space for endogenous variable p
Np              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
pgrid           = nodeunif(Np,glob.pmin.^curv(1),glob.pmax.^curv(1)).^(1/curv(1));  % Adds curvature
pgrid0          = pgrid;    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',pgrid,0,glob.spliorder(1)},...
                         {'spli',agrid,0,glob.spliorder(2)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
pgrid           = s(s(:,2)==s(1,2),1); 
agrid           = s(s(:,1)==s(1,1),2);
Np              = size(pgrid,1); 
Na              = size(agrid,1);

%% Compute expectations matrix
Phi             = funbas(fspace,s);
glob.Emat       = kron(A,speye(Np))*Phi;

%% Construct fine grid for histogram
pgridf          = nodeunif(glob.nf(1),glob.pmin.^glob.curv(1),glob.pmax.^glob.curv(1)).^(1/glob.curv(1));
Npf             = size(pgridf,1);

agridf          = agrid;     
Naf             = size(agridf,1);
sf              = gridmake(pgridf,agridf);
Nsf             = size(sf,1);

glob.pgridf     = pgridf;
glob.agridf     = agridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QA matrix for approximation of stationary distribution
glob.QA         = kron(A,ones(Npf,1)); 

%% Create one time only basis matrices: these might need to be changed
glob.Phi_A      = splibas(agrid0,0,spliorder(2),s(:,2));                % Used in Bellman / Newton computing expected values
glob.Phi_Af     = splibas(agrid0,0,spliorder(2),sf(:,2));               % Used when solving on fine grid
Phi_P           = splibas(pgrid0,0,spliorder(1),s(:,1)*1/exp(param.mu));
glob.Phi        = dprod(glob.Phi_A,Phi_P);                              % Used in Bellman / Newton updating of c

%% Declare additional global variables
glob.pgrid0     = pgrid0;
glob.pgrid      = pgrid;
glob.agrid0     = agrid0;
glob.agrid      = agrid;
glob.A          = A;
glob.Assa       = Assa;
glob.Na         = Na;
glob.Np         = Np;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;

end