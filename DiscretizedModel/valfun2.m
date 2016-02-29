function [ valout ] = valfun2(flag,parms,pP,V0,prof)

%valfun2(flag,parms,pP,V0,prof)
%VALFUN Summary of this function goes here
%-------------------------------------------------
%   Compute the firm's value function as function of past value function,
%   parameters, the old relative price, productivity, 
%   aggregate output/demand, money, inflation, and decision about whether 
%   to pay menu cost.
%
%   Do this one state-space element at a time, because different points in
%   the space lead to either 'K' or 'C'
%
%   INPUTS
%   - flag  = 'K' if keep old price, 'C' if change to new price
%   - parms = structure; model parameters
%       * .Ppgrid = vector; grid points for real prices (Npp x 1) 
%       * .grid  = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i   = scalar; size of particular state variable vector e.g. Na
%       * .trans = matrix; transition prob matrix (Na*Ndm*Ny*Npi x Na*Ndm*Ny*Npi)
%   - pP    = current real price
%   - V0    = matrix; previous iteration's value function (Na*Ndm*Ny*Npi x Np)
%
%-------------------------------------------------

beta  = parms.beta;
pPmin = parms.pPmin;
pPmax = parms.pPmax;

Npp   = parms.Npp;
Na    = parms.Na;
Ny    = parms.Ny;
Ndm   = parms.Ndm;
% gridsize = Na*Ny*Ndm;

grid = parms.grid;
pPgrid = parms.pPgrid;
trans = parms.trans;

% Compute inflation rate tomorrow given every state variable today, and 
% every state variable tomorrow 
% in logs: pi_t+1 = Dm_t+1 + ln Y_t+1 - ln Y_t
% 


switch flag
    case 'K'  % Keep old price
        E_V = nan(Na*Ndm*Ny,1);
        
        for i = 1:Na*Ndm*Ny
            PI =  (grid(3,:).*grid(2,:))/grid(2,i);
            tmp = abs( repmat(pPgrid,Na*Ndm*Ny,1) - repmat( pP./PI',1,Npp) ); % how far from grid?
            [~, idx_pP] = min(tmp,[],2); % index of closest values
            E_V(i,1) = V0(i,idx_pP(i));  % Move to new pP grid point in V0
        end
        
%         E_V = nan(Na*Ndm*Ny,1);       
%         PI =  grid(2,:)'*(grid(3,:).*grid(2,:)).^(-1); % pi_t+1 = Dm_t+1 + ln Y_t+1 - ln Y_t   
%         ptmp = reshape(pPgrid,1,1,Npp);
%         tmp = abs( repmat(ptmp,Na*Ndm*Ny,Na*Ndm*Ny,1) - repmat(pP./PI',1,1,Npp) ); % how far from grid?
%         [~, idx_pP] = min(tmp,[],3); % index of closest values
% 
%         for i = 1:Na*Ndm*Ny
%             for j = 1:Na*Ndm*Ny
%             % Find new value for tomorrow's price: pP/(P'/P) = pP'            
%             %  PI =  (grid(3,:).*grid(2,:))/grid(2,i);
% 
%             E_V(i,j) = V0(i,idx_pP(i,j));  % Move to new pP grid point in V0                       
%             end
%         end
        valout = prof + beta*trans*E_V;
        
    case 'C' % Change to new price
        
        E_V = nan(Na*Ndm*Ny,Npp);
        for i = 1:Na*Ndm*Ny
            PI =  (grid(3,:).*grid(2,:))/grid(2,i);
            for j = 1:Npp               
                tmp = abs( repmat(pPgrid,Na*Ndm*Ny,1) - repmat( pPgrid(1,j)./PI',1,Npp) ); % how far from grid?
                [~, idx_pP] = min(tmp,[],2); % index of closest values                
                E_V(i,j) = V0(i,idx_pP(i));  % Move to new pP grid point in V0
            end
        end
        val2max =  prof  + beta*trans*E_V;
%         val2max =  realprofit('C',parms,pPgrid)  + beta*trans*E_V;

        [valout, ~] = max(val2max,[],2); % max value and pPgrid index of max value
end




end

