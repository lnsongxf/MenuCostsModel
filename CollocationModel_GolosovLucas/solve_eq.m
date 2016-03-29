function [eq] = solve_eq(param,glob,options)
%SOLVE_EQ Solve for equilibrium in the no aggregate uncertainty case
%-------------------------------------------------
%   Solves for Y and P that satisfy equilibrium when there is no aggregate
%   uncertainty in the model.
%
%   INPUTS
%   - Y         = conjectured value of output, Y
%   - P         = conjectured value of aggregate price level, P
%   - param     = parameters 
%   - glob      = includes state space, function space, approximating functions etc
%   - options   = 
%   OUTPUT
%   - eq        = equilibrium objects: Y, P
%-------------------------------------------------

% initialize
Yout = 1/2*(options.Ylb + options.Yub);

% Storage
Yinvec          = zeros(options.itermaxp,1); 
Youtvec         = zeros(options.itermaxp,1);
options.cresult = [];
tictic = tic;

for tt = (1:options.itermaxp)
    % 1. Update Y 
    Yin             = Yout;
    % 2. Solve economy given Y
    eq              = solve_cL(Yin,param,glob,options);  
    options.cresult = eq.c;         % Save to use as starting guess
    Yout            = glob.damp*Yin   + (1-glob.damp)*eq.Y;     % update Y
    % 3. Record output and print
    Yinvec(tt)      = Yin;
    Youtvec(tt)     = eq.Y;
    if strcmp(options.eqprint,'Y') 
        fprintf('%2i. Yin:\t%2.4f\tYout:\t%2.4f\tt:%2.1f\n',tt,Yin,eq.Y,toc(tictic));
    end
    % 6. Break if equilibrium
    if abs(Yin-eq.Y)<options.tolYeq
        fprintf('Solved for output, t = %1i\n', tt)
        break
    end  
    % 7. Plot option
    if strcmp(options.eqplot,'Y') 
       figure(options.fignum);
       subplot(2,2,3);
       plot(Yinvec(Yinvec~=0),'bo','markersize',6,'markerfacecolor','b');hold on;grid on;
       plot(Youtvec(Youtvec~=0),'ro','markersize',6,'markerfacecolor','r');
       xlim([0,max(tt,12)]);
       legend('Yin','Yout','Location','NorthEast');
       title('Equilibrium');
       set(gca,'fontsize',options.fontsize)
       drawnow;  
    end
end
