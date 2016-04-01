function [param,glob] = setup_ss(param,glob,options)

%% State space for idiosyncratic productivity
% One persistent shock
Nnu                 = glob.n(2); 
[Nu,nugrid,Nussa]   = setup_MarkovZ(Nnu,param.sigmanu,(1-param.eta),1);         
nugrid0             = nugrid;

%% State space for endogenous variable x
Nx              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
xgrid           = nodeunif(Nx,glob.xmin.^curv(1),glob.xmax.^curv(1)).^(1/curv(1));  % Adds curvature
xgrid0          = xgrid;    % Save for computing basis matrices

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',xgrid,0,glob.spliorder(1)},...
                         {'spli',nugrid,0,glob.spliorder(2)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
xgrid           = s(s(:,2)==s(1,2),1); 
nugrid          = s(s(:,1)==s(1,1),2);
Nx              = size(xgrid,1); 
Nnu             = size(nugrid,1);

%% Compute expectations matrix
Phi             = funbas(fspace,s);
glob.Emat       = kron(Nu,speye(Nx))*Phi;

%% Construct fine grid for histogram
xgridf          = nodeunif(glob.nf(1),glob.xmin.^glob.curv(1),glob.xmax.^glob.curv(1)).^(1/glob.curv(1));
Nxf             = size(xgridf,1);

nugridf         = nugrid;     
Nnuf            = size(nugridf,1);
sf              = gridmake(xgridf,nugridf);
Nsf             = size(sf,1);

glob.xgridf     = xgridf;
glob.nugridf    = nugridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QNu matrix for approximation of stationary distribution
glob.QNu         = kron(Nu,ones(Nxf,1)); 

%% Create one time only basis matrices
glob.Phi_nu     = splibas(nugrid0,0,spliorder(2),s(:,2));                % Used in Bellman / Newton computing expected values
glob.Phi_nuf    = splibas(nugrid0,0,spliorder(2),sf(:,2));               % Used when solving on fine grid
% check this!!!!!!
Phi_X           = splibas(xgrid0,0,spliorder(1),s(:,1)*1/exp(param.mu));
glob.Phi        = dprod(glob.Phi_nu,Phi_X);                              % Used in Bellman / Newton updating of c

%% Declare additional global variables
glob.xgrid0     = xgrid0;
glob.xgrid      = xgrid;
glob.nugrid0    = nugrid0;
glob.nugrid     = nugrid;
glob.Nu         = Nu;
glob.Nussa      = Nussa;
glob.Nnu        = Nnu;
glob.Nx         = Nx;
glob.fspace     = fspace;
glob.s          = s;
glob.Ns         = Ns;

end