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

Ylb0    = options.Ylb;
Yub0    = options.Yub;
Ylb     = Ylb0;
Yub     = Yub0;

% Storage
Yinvec          = zeros(options.itermaxp,1); 
Youtvec         = zeros(options.itermaxp,1);
options.cresult = [];
tictic = tic;
for tt = (1:options.itermaxp)
    % 1. Update Y 
    Y               = (1/2)*(Ylb+Yub);
    % 2. Solve economy given Y
    eq              = solve_cL(Y,param,glob,options);  
    options.cresult = eq.c;         % Save to use as starting guess
    % 3. Record output and print
    Yinvec(tt)      = Y;
    Youtvec(tt)     = eq.Y;
    if strcmp(options.eqprint,'Y') 
        fprintf('%2i. Yin:\t%2.4f\tYout:\t%2.4f\tt:%2.1f\n',tt,Y,eq.Y,toc(tictic));
%         fprintf('%2i. Pin:\t%2.4f\tPout:\t%2.4f\tt:%2.1f\n',tt,1,eq.P,toc(tictic));
    end
    % 4. Set all flags
    d               = Yinvec-Youtvec;
    eq.flag.exist   = ~all(sign(d)==max(sign(d))); 
    eq.flag.equi    = (abs(Yinvec(tt)-Youtvec(tt))<options.tolYeq);
    eq.flag.down    = (Yinvec(tt)>Youtvec(tt));
    eq.flag.up      = (Yinvec(tt)<Youtvec(tt));
    % 5. Shift bounds
    Ylb             = (eq.flag.up)*Y    + (eq.flag.down)*Ylb;
    Yub             = (eq.flag.up)*Yub  + (eq.flag.down)*Y;
    % 6. Break if equilibrium
    if eq.flag.equi
        break
    end
    % 7. Plot option
    if strcmp(options.eqplot,'Y') 
       figure(options.fignum)       
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
