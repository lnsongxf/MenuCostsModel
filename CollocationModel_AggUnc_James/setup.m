function [param,glob] = setup_agg(param,glob,cKS,options)
%SETUP Prepares model objects for collocation procedure
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

switch options.AR1
%     case 'Y' 
%         logzub              = norminv(1-glob.pzlb,0,glob.sige/sqrt(1-glob.rhoz^2));
%         logzlb              = norminv(glob.pzlb,0,glob.sige/sqrt(1-glob.rhoz^2));
%         zlb                 = exp(logzlb);
%         zub                 = exp(logzub);
%         agrid               = nodeunif(Na,zlb,zub);
%         % For comparison compute Rouwenhurst grid and Pssz
%         [P,agridRouw,Pssz]  = setup_MarkovZ(Na,glob.sige,glob.rhoz,1);
    case 'N'
        [Pa,agrid,Pssa]      = setup_MarkovZ(Na,param.sigzeta,param.rhoa,1);
        
        muvar = [param.mu*(1-param.rhom); b0 + b2*param.mu*(1-param.rhom)];
        Avar = [param.rhom, 0; b2*param.rhom, b1];
        Svar = [param.sigmaeps; b2*param.sigmaeps];
        [tmpgrid,Pmy,Pssmy]              = tauchenvar([Nm; Ny],muvar,Avar,Svar);
        Pmy = Pmy';     % Make sure 
        Mgrid = exp(unique(tmpgrid(1,:)))'; % Take D(ln(M_t)) back to DM_t 
        Ygrid = exp(unique(tmpgrid(2,:)))'; % Take ln(Y_t) back to Y_t         
        
end
agrid0 = agrid;
Mgrid0 = Mgrid;
Ygrid0 = Ygrid;
%% State space for endogenous variable pP
NpP              = glob.n(1);
curv            = glob.curv;
spliorder       = glob.spliorder;
pPgrid           = nodeunif(NpP,glob.pPmin.^curv(1),glob.pPmax.^curv(1)).^(1/curv(1));  % Adds curvature
pPgrid0          = pPgrid;    % Save for computing basis matrices in valfunc.m:line9

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',pPgrid,0,spliorder(1)},...
                         {'spli',agrid,0,spliorder(2)},...
                         {'spli',Mgrid,0,spliorder(3)},...
                         {'spli',Ygrid,0,spliorder(4)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
pPgrid = sgrid{1};  %s(s(:,2)==s(1,2),1); 
agrid = sgrid{2};   % s(s(:,1)==s(1,1),2);
Mgrid = sgrid{3};
Ygrid = sgrid{4};

NpP = size(pPgrid,1); 
Na = size(agrid,1);
Nm = size(Mgrid,1);
Ny = size(Ygrid,1);

%% Compute expectations matrix
switch options.AR1
%     case 'Y'
%         Ne              = glob.Ne1;
%         pvec            = nodeunif(Ne,glob.plb,1-glob.plb);     % Make an equi-spaced grid in probabilities
%         e               = norminv(pvec,0,glob.sige);            
%         w               = normpdf(e,0,glob.sige);               % Invert normal for shocks
%         w               = w/sum(w);                             % Compute pdf of shocks
%         iNe             = ones(Ne,1);                           
%         iNs             = ones(Ns,1);
%         gfun            = @(z,e) max(min(exp(glob.rhoz*log(z)+e),max(agrid)),min(agrid));   % Constrained to lie within nodes
%         g               = gfun(kron(s(:,2),iNe),kron(iNs,e));
%         Phi             = funbas(fspace,[kron(s(:,1),iNe),g]);
%         Ikronw          = kron(eye(Ns),w');
%         glob.Emat       = Ikronw*Phi;        
    case 'N'
        Phi             = funbas(fspace,s);
        
        
% NOTE: This needs to be post-(dot)multiplied by an indicator function before being
% multilplied by a collaction coefficient vector. 


        H          = kron(ones(1,Ns),kron(speye(Ny),ones(NpP*Na*Nm,Ny*Nm)));
        glob.H          = bsxfun(@rdivide,H,sum(H,2));    % Normalize rows of H to sum to one
        glob.Emat  = kron(kron(Pmy,Pa),speye(NpP));
%         glob.Emat       = kron(Pa,speye(NpP));      % In no agg uncertainty case   
        
end

%% Construct fine grid for histogram
pPgridf          = nodeunif(glob.nf(1),glob.pPmin.^glob.curv(1),glob.pPmax.^glob.curv(1)).^(1/glob.curv(1));
NpPf             = size(pPgridf,1);

switch options.AR1
%     case 'Y'
%         agridf  = nodeunif(glob.nf(2),min(agrid),max(agrid));
    case 'N'
        agridf  = agrid;
end        

Naf             = size(agridf,1);
sf              = gridmake(pPgridf,agridf,Mgrid,Ygrid);
Nsf             = size(sf,1);

glob.pPgridf     = pPgridf;
glob.agridf     = agridf;
glob.sf         = sf;
glob.Nsf        = Nsf;

%% Compute QZ matrix for approximation of stationary distribution
switch options.AR1
%     case 'Y' 
%         % 1. Continuous AR(1)
%         Ne              = glob.Ne2;
%         pvec            = nodeunif(Ne,glob.plb,1-glob.plb);         % Make an equi-spaced grid in probabilities
%         e               = norminv(pvec,0,glob.sige);                % Invert normal for shocks
%         w               = normpdf(e,0,glob.sige);                   % Compute pdf of shocks
%         w               = w/sum(w);                                 % Normalise
%         fspaceZ         = fundef({'spli',agridf,0,1});              % Linear interpolant
%         QZ              = zeros(Nsf,Naf);
%         P               = zeros(Naf,Naf);                           % P constructed so can compute steady state Psszf and compare to Pssz
%         for i = 1:Ne;
%             g           = gfun(sf(:,2),e(i));
%             QZi         = funbas(fspaceZ,g);
%             QZ          = QZ + w(i)*QZi;
%             P           = P  + w(i)*funbas(fspaceZ,gfun(agridf,e(i)));
%         end
%         glob.QZ         = QZ;
%         % For plotting
%         Psszf           = P^1000;
%         Psszf           = Psszf(1,:)';
%         glob.Psszf      = Psszf;
%         if strcmp(options.plotSD,'Y')
%             figure(round(1000*rand));
%             subplot(1,2,1);
%             plot(agridRouw,Pssz,'bo-');grid on; hold on;title('A. Coarse Transition matrix Pssz - From Rouwenhurst'); 
%             subplot(1,2,2);
%             plot(agridf,Psszf,'ro-');grid on;title('B. Fine stationary dist');
%         end
    case 'N'
        % 2. Discrete productivity process
        glob.QA         = kron(Pa,ones(NpPf,1)); 
end

%% Create basis matrices (need only compute once)
glob.Phi_A      = splibas(agrid0,0,spliorder(2),s(:,2));        % Used in Bellman / Newton computing expected values
glob.Phi_Af     = splibas(agrid0,0,spliorder(2),sf(:,2));       % Used when solving on fine grid
Phi_pP          = splibas(pPgrid0,0,spliorder(1),s(:,1));       % Basis matrix for price grid
Phi_pPf         = splibas(pPgrid0,0,spliorder(1),sf(:,1));      % Used when solving on fine grid
glob.Phi_M      = splibas(Mgrid0,0,spliorder(3),s(:,3));        % Used in Bellman / Newton computing expected values
glob.Phi_Y      = splibas(Ygrid0,0,spliorder(4),s(:,4));        % Used in Bellman / Newton computing expected values
glob.Phi_Mf      = splibas(Mgrid0,0,spliorder(3),sf(:,3));        % Used in Bellman / Newton computing expected values
glob.Phi_Yf      = splibas(Ygrid0,0,spliorder(4),sf(:,4));        % Used in Bellman / Newton computing expected values

% [Np,Na,Nm,Ny] --> Work backwards inside the dprod function
glob.Phi        = dprod(glob.Phi_Y,dprod(glob.Phi_M,dprod(glob.Phi_A,Phi_pP)));                     % Used in Bellman / Newton updating of c

% Don't need to put aggregate states on a fine grid.
glob.Phif        = dprod(glob.Phi_Yf,dprod(glob.Phi_Mf,dprod(glob.Phi_Af,Phi_pPf)));                  % Used when solving on fine grid

%% Create complicated basis matrix: 
% Phi(s1 kron Ygrid.^-1 kron Ygrid' kron Mgrid'.^-1, s kron i_{NyNyNm) )

s_prime = [kron(kron(kron(s(:,1), Ygrid.^(-1)), Ygrid), Mgrid.^(-1)), kron(s(:,2:end), ones(Ny*Ny*Nm,1))];

Phi_pPprime  = splibas(pPgrid0,0,spliorder(1),s_prime(:,1));       % Basis matrix for price grid
Phi_Aprime   = splibas(agrid0,0,spliorder(2),s_prime(:,2));        % Used in Bellman / Newton computing expected values
Phi_Mprime   = splibas(Mgrid0,0,spliorder(3),s_prime(:,3));        % Used in Bellman / Newton computing expected values
Phi_Yprime   = splibas(Ygrid0,0,spliorder(4),s_prime(:,4));        % Used in Bellman / Newton computing expected values

% [Np,Na,Nm,Ny] --> Work backwards inside the dprod function
glob.Phiprime        = dprod(Phi_Yprime, dprod(Phi_Mprime, dprod(Phi_Aprime, Phi_pPprime)));       


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
glob.Pssa       = Pssmy;             % stationary distribution for M,Y
glob.Na         = Na;               % length of state space grid for a 
glob.NpP        = NpP;              % length of state space grid for pP
glob.Nm         = Nm;               % length of state space grid for M 
glob.Ny         = Ny;               % length of state space grid for Y 
glob.fspace     = fspace;           % function space object for the model
glob.s          = s;                % full state space 
glob.Ns         = Ns;               % size of full state space
glob.s_prime    = s_prime;          % Adjusted price 


    
    
    
    
end

        
   
