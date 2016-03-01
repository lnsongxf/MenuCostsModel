function [V_final, Vc_final, Vk_final] = val_iter( parms, tol)
%VAL_ITER Implement value function iteration on the menu costs model
%-------------------------------------------------
%   Implement value function iteration on the menu costs model.
%
%   INPUTS
%   - parms      = structure; model parameters
%       * .pPgrid= matrix; grid points for relative price [Npp]
%       * .grid  = matrix; grid points for model: [Na;Ndm;Ny;Npi];
%       * .N_i   = scalar; size of particular state variable vector e.g. Na
%   - T          = Number of iterations
%
%   OUTPUTS
%   - V_final   = The solved value function V = max{Vk,Vc}
%   - Vc_final  = The solved value function if changing price
%   - Vk_final  = The solved value function if keeping current price
%--------------------------------------------------
% NOTES: MUCH faster now. 
% 30,000 grid points ~ 1 minute
% 60,000 grid points ~ 2 minutes.
% 160,000 grid points ~ 9 minutes

% Can make even faster by eliminating the 'j' loop.

% Precompute expectations for additional speed? Check accuracy, seems to be
% something wrong when we do that...

% Look for other speed improvements...
%--------------------------------------------------------------

% Initialize Value function matrices
Vk(:,:) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);
Vc(:,:) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);
V(:,:) = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);

% Load previous solutionss
if exist(['FinalV_Na' num2str(parms.Na) '_Ndm' num2str(parms.Ndm) '_Ny' num2str(parms.Ny) '_Npi' num2str(parms.Npi) '_Npp' num2str(parms.Npp)], 'file') == 2   
    tmp = load(['FinalV_Na' num2str(parms.Na) ...
        '_Ndm' num2str(parms.Ndm) ...
        '_Ny' num2str(parms.Ny) ...
        '_Npi' num2str(parms.Npi) ...
        '_Npp' num2str(parms.Npp)]);
    V(:,:) = tmp.V_final;
else  
    V(:,:) = eye(parms.Na*parms.Ny*parms.Ndm*parms.Npi, parms.Npp);
end


% pre-allocate real profit functions for speed
profK = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi, parms.Npp);
for j = 1:parms.Npp
    % realprofit(flag,parms,pP,a,Y)
    profK(:,j)  = ...
        realprofit('K',parms,parms.pPgrid(j),parms.grid(1,:),parms.grid(3,:));
end
profC = realprofit('C',parms,parms.pPgrid,parms.grid(1,:),parms.grid(3,:));


% % Pre-find matrix indicies when changing
% for j = 1:parms.Npp
%     tmp = abs( repmat(parms.pPgrid,parms.Na*parms.Ndm*parms.Ny*parms.Npi,1) - ...
%         repmat(parms.pPgrid(1,j)./parms.grid(4, :)',1,parms.Npp) ); % how far from grid?
%     [~, idx_pP] = min(tmp,[],2); % index of closest values
% end

V_prev = ones(parms.Na*parms.Ny*parms.Ndm*parms.Npi, parms.Npp);
while norm(V - V_prev) > tol
    
    % Precompute expectation for speed
%     E_V = nan(parms.Na*parms.Ndm*parms.Ny*parms.Npi,parms.Npp);
%     for i = 1:parms.Na*parms.Ndm*parms.Ny*parms.Npi
%         E_V(i,j) = V(i,idx_pP(i),t);  % Move to new pP grid point in V0
%     end
    
    Vc(:,:)   = repmat(valfun('C', parms, parms.pPgrid, V(:,:),profC), 1, parms.Npp);
    
    for j = 1:parms.Npp 
        Vk(:,j)   = valfun('K', parms, parms.pPgrid(j), V(:,:), profK(:,j));        
    end
    
    V_prev = V;
    V(:,:) = bsxfun(@max,Vk(:,:),Vc(:,:));  % bsxfun = binary operations on matrices

    disp(['Norm = ' num2str(norm(V - V_prev)) ]) 
end

V_final = V(:,:);
Vc_final = Vc(:,:);
Vk_final = Vk(:,:);

% Save the final value fuction 
save(['FinalV_Na' num2str(parms.Na) ...
    '_Ny' num2str(parms.Ny) ... 
    '_Ndm' num2str(parms.Ndm) ...
    '_Npi' num2str(parms.Npi) ...
    '_Npp' num2str(parms.Npp) '.mat'], 'V_final');






end

