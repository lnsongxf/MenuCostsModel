function [ pPout ] = newprice(parms,pP,statenum,Vk,Vc,V)
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
%   - flag          = 'K' if keep old price, 'C' if change to new price
%   - parms         = structure; model parameters
%       * .Ppgrid   = vector; grid points for real prices (Npp x 1) 
%       * .grid     = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i      = scalar; size of particular state variable vector e.g. Na
%       * .trans    = matrix; transition prob matrix (Na*Ndm*Ny*Npi x Na*Ndm*Ny*Npi)
%   - pP            = current real price
%   - statenum      = index number representing the current state 
%   - V             = matrix; value function (Na*Ndm*Ny*Npi x Np)
%
%   OUTPUTS
%   - pPout = new price. If 'K', same as before. If 'C', new argmax.
%-------------------------------------------------

a = parms.grid(1,statenum);
Y = parms.grid(3,statenum);
trans = parms.trans(statenum,:);
pPgrid = parms.pPgrid;

% Check whether firm keeps or changes price
tmp = abs( pPgrid - pP );   % how far from grid?
[~, idx_pP] = min(tmp);     % index of closest values   
pP = pPgrid(idx_pP);       % Put price on the grid

if Vk(statenum,indx_pP) > Vc(statenum,indx_pP)

    %      flag = 'K';
    pPout = pP;

elseif Vk(statenum,indx_pP) <= Vc(statenum,indx_pP)

    %     flag = 'C';
    val2max =  realprofit('C',parms,pP,a,Y)  + beta*trans*V; 
    [~, idx_new] = max(val2max); % max value and pPgrid index of max value
    pPout = pPgrid(idx_new,1);
    
end




end


