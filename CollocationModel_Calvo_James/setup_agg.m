function [glob] = setup_agg(param,glob,cKS,options)
%SETUP_AGG Prepares model objects for collocation procedure, agg uncertainty case
%-------------------------------------------------
%   This file prepares the Markov process (continuous or discretized), the
%   state space (across exogenous and endogenous variables), the function
%   space for approximating with, elements for use in approximation of the
%   stationary distribution, and basis functions for use in the collocation
%   procedure computations
%
%   INPUTS
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - cKS       = Krussel Smith law of motion parameters: (lnY = b0 + b1lnY-1 + b2 Dm)
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%-------------------------------------------------
%% State space for idiosyncratic productivity
% One persistent shock  % [Np,Na,Nm,Ny]
Na              = glob.n(2);
Nm              = glob.n(3);
Ny              = glob.n(4);
b0 = cKS(1);
b1 = cKS(2);
b2 = cKS(3);

[Pa,agrid,Pssa]      = setup_MarkovZ(Na,param.sigzeta,param.rhoa,1);

A0 = eye(2,2);
A1 = [param.mu*(1-param.rhom); b0 + b2*param.mu*(1-param.rhom)];
A2 = [param.rhom, 0; b2*param.rhom, b1];
SIGMA = [param.sigmaeps^2, param.sigmaeps^2*b2;...
          param.sigmaeps^2*b2, param.sigmaeps^2*b2^2];
N = [Nm; Ny];
method = 1;    % Method=1: uniform grid, Method=2: uses normal CDF to pick grid

[Pmy,my_grid,~] = fn_var_to_markov(A0,A1,A2,SIGMA,N,1000,method);
% Test the Markov chain
% EZ = (eye(2,2) - A2)^(-1)*A1;    % Theoretical mean of VAR
% VarZ = lyap(A2,SIGMA);           % Theoretical variance
Mgrid = exp(unique(my_grid(1,:)))'; % Take D(ln(M_t)) back to DM_t
Ygrid = exp(unique(my_grid(2,:)))'; % Take ln(Y_t) back to Y_t
        
agrid0 = agrid;
Mgrid0 = Mgrid;
Ygrid0 = Ygrid;
%% State space for endogenous variable pP
NpP             = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
pPgrid          = nodeunif(NpP,glob.pPmin.^curv(1),glob.pPmax.^curv(1)).^(1/curv(1));  % Adds curvature
pPgrid0         = pPgrid;    % Save for computing basis matrices in valfunc.m:line9

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',pPgrid,0,spliorder(1)},...
                         {'spli',agrid,0,spliorder(2)},...
                         {'spli',Mgrid,0,spliorder(3)},...
                         {'spli',Ygrid,0,spliorder(4)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
pPgrid          = sgrid{1};  
agrid           = sgrid{2};   
Mgrid           = sgrid{3};
Ygrid           = sgrid{4};

NpP             = size(pPgrid,1); 
Na              = size(agrid,1);
Nm              = size(Mgrid,1);
Ny              = size(Ygrid,1);

%% Compute expectations matrix
glob.Phi        = funbas(fspace,s);
H               = kron(ones(1,Ns),kron(speye(Ny),ones(NpP*Na*Nm,Ny*Nm)));
glob.H          = bsxfun(@rdivide,H,sum(H,2));    % Normalize rows of H to sum to one
glob.Emat       = kron(kron(Pmy,Pa),speye(NpP));


%% Construct fine grid for histogram
pPgridf         = nodeunif(glob.nf(1),glob.pPmin.^glob.curv(1),glob.pPmax.^glob.curv(1)).^(1/glob.curv(1));
NpPf            = size(pPgridf,1);
agridf          = agrid;
Naf             = size(agridf,1);
sf              = gridmake(pPgridf,agridf,Mgrid,Ygrid);
Nsf             = size(sf,1);

glob.pPgridf    = pPgridf;
glob.agridf     = agridf;
glob.sf         = sf;
glob.Nsf        = Nsf;
glob.Phif       =funbas(fspace,sf);

%% Compute QA matrix for approximation of stationary distribution
glob.QA         = kron(Pa,ones(NpPf,1)); 

%% Create basis matrices (need only compute once)
glob.Phi_A      = splibas(agrid0,0,spliorder(2),s(:,2));        % Used in Bellman / Newton computing expected values
glob.Phi_Af     = splibas(agrid0,0,spliorder(2),sf(:,2));       % Used when solving on fine grid
glob.Phi_M      = splibas(Mgrid0,0,spliorder(3),s(:,3));        % Used in Bellman / Newton computing expected values
glob.Phi_Mf     = splibas(Mgrid0,0,spliorder(3),sf(:,3));        % Used in Bellman / Newton computing expected values
glob.Phi_Y      = splibas(Ygrid0,0,spliorder(4),s(:,4));        % Used in Bellman / Newton computing expected values
glob.Phi_Yf     = splibas(Ygrid0,0,spliorder(4),sf(:,4));        % Used in Bellman / Newton computing expected values

%% Create the more complicated basis matrix: 
% Phi(s1 kron Ygrid.^-1 kron Ygrid' kron Mgrid'.^-1, s kron i_{NyNyNm) )
s_prime         = [kron(kron(kron(s(:,1), Ygrid.^(-1)), Ygrid), Mgrid.^(-1)), kron(s(:,2:end), ones(Ny*Ny*Nm,1))];
glob.Phiprime   = funbas(fspace,s_prime);

%% Declare additional global variables
glob.pPgrid0    = pPgrid0;          % unique elements of pP grid
glob.pPgrid     = pPgrid;           % full state space for pP grid
glob.agrid0     = agrid0;           % unique elements of a grid
glob.agrid      = agrid;            % full state space for a grid
glob.Mgrid0     = Mgrid0;           % unique elements of a grid
glob.Mgrid      = Mgrid;            % full state space for a grid
glob.Ygrid0     = Ygrid0;           % unique elements of a grid
glob.Ygrid      = Ygrid;            % full state space for a grid
glob.Pa         = Pa;              % transition probability matrix for a
glob.Pssa       = Pssa;             % stationary distribution for a
glob.Pmy        = Pmy;              % transition probability matrix for M,Y
% glob.Pssmy      = Pssmy;            % stationary distribution for M,Y
glob.Na         = Na;               % length of state space grid for a 
glob.NpP        = NpP;              % length of state space grid for pP
glob.Nm         = Nm;               % length of state space grid for M 
glob.Ny         = Ny;               % length of state space grid for Y 
glob.fspace     = fspace;           % function space object for the model
glob.s          = s;                % full state space 
glob.Ns         = Ns;               % size of full state space
glob.s_prime    = s_prime;          % Adjusted price 
    
end

        
   
