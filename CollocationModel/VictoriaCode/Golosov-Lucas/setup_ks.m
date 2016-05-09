function [param,glob] = setup_ks(cKS,param,glob,options)

%% State space for idiosyncratic productivity
% One persistent shock
Nnu                 = glob.n(2); 
[Nu,nugrid,Nussa]   = setup_MarkovZ(Nnu,param.sigmanu,(1-param.eta),1);         
nugrid0             = nugrid;

%% State space for aggregate cbar
Nc                  = glob.n(3);
c_mean              = (cKS.b0 + cKS.b2*param.mu)/(1-cKS.b1);
[C,cgrid,Cssa]      = setup_MarkovZ(Nc,param.sigmam*cKS.b2,cKS.b1,exp(c_mean));
cgrid0              = cgrid;

%% State space for endogenous variable x
Nx              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
xgrid           = nodeunif(Nx,glob.xmin.^curv(1),glob.xmax.^curv(1)).^(1/curv(1));  % Adds curvature
xgrid0          = xgrid;    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',xgrid,0,glob.spliorder(1)},...
                         {'spli',nugrid,0,glob.spliorder(2)},...
                         {'spli',cgrid,0,glob.spliorder(3)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
%xgrid           = s(s(:,2)==s(1,2),1); 
%nugrid          = s(s(:,1)==s(1,1),2);
xgrid           = unique(s(:,1));
nugrid          = unique(s(:,2));
cgrid           = unique(s(:,3));
Nx              = size(xgrid,1); 
Nnu             = size(nugrid,1);
Nc              = size(cgrid,1);

%% Discrete pdf and support for money shocks (quadrature)
[glob.supp_e glob.f]      = qnwnorm(glob.Ne,0,param.sigmam^2);

%% Compute expectations matrix
%Phi             = funbas(fspace,s);
glob.Emat       = kron(kron(C,Nu),speye(Nx))*kron(speye(Nx*Nnu*Nc),glob.f');

%% Construct fine grid for histogram
xgridf          = nodeunif(glob.nf(1),glob.xmin.^glob.curv(1),glob.xmax.^glob.curv(1)).^(1/glob.curv(1));
Nxf             = size(xgridf,1);

nugridf         = nugrid;     
Nnuf            = size(nugridf,1);

cgridf          = cgrid;     
Ncf             = size(cgridf,1);

sf              = gridmake(xgridf,nugridf,cgridf);
Nsf             = size(sf,1);

glob.xgridf     = xgridf;
glob.nugridf    = nugridf;
glob.cgridf     = cgridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QNu matrix for approximation of stationary distribution
glob.QNu       = kron(Nu,ones(Nxf,1)); 

%% Create one time only basis matrices
glob.Phi_nu     = splibas(nugrid0,0,spliorder(2),s(:,2));                % Used in Bellman / Newton computing expected values
glob.Phi_nuf    = splibas(nugrid0,0,spliorder(2),sf(:,2));               % Used when solving on fine grid

glob.Phi_c      = splibas(cgrid0,0,spliorder(3),s(:,3));                % Used in Bellman / Newton computing expected values
glob.Phi_cf     = splibas(cgrid0,0,spliorder(3),sf(:,3));               % Used when solving on fine grid

% construct continuation value
x_col           = kron(s(:,1),ones(glob.Ne,1)).*1/exp(param.mu).*(1./kron(ones(Nx*Nnu*Nc,1),exp(glob.supp_e)));
nuc_col         = kron(s(:,2:end),ones(glob.Ne,1));
glob.Phi_stilde = funbas(fspace,[x_col nuc_col]);

Phi_X           = splibas(xgrid0,0,spliorder(1),s(:,1));
glob.Phi        = dprod(glob.Phi_c,dprod(glob.Phi_nu,Phi_X));                              % Used in Bellman / Newton updating of c

%% Declare additional global variables
glob.xgrid0     = xgrid0;
glob.xgrid      = xgrid;
glob.nugrid0    = nugrid0;
glob.nugrid     = nugrid;
glob.cgrid0     = cgrid0;
glob.cgrid      = cgrid;
glob.Nu         = Nu;
glob.Nussa      = Nussa;
glob.Nnu        = Nnu;
glob.C          = C;
glob.Cssa       = Cssa;
glob.Nc         = Nc;
glob.Nx         = Nx;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;

end