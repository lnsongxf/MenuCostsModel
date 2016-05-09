function [eq,val] = solve_eq_menucost(param,glob,options)

Ylb0    = options.Ylb;
Yub0    = options.Yub;
Ylb     = Ylb0;
Yub     = Yub0;

% Storage
Yinvec          = zeros(options.itermaxY,1); 
Youtvec         = zeros(options.itermaxY,1);
options.cresult = [];
tictic = tic;

for tt = (1:options.itermaxY)
    % 1. Update p 
    Y               = (1/2)*(Ylb+Yub);
    % 2. Solve economy given p
    eq              = solve_pL(Y,param,glob,options);  
    options.cresult = eq.c;         % Save to use as starting guess
    % 3. Record output and print
    Yinvec(tt)      = Y;
    Youtvec(tt)     = eq.G_Y;
    if strcmp(options.eqprint,'Y') 
        fprintf('%2i. Yin:\t%2.6f\tdistance:\t%2.6f\tt:%2.1f\n',tt,Y,eq.G_Y,toc(tictic));
    end
    % 4. Set all flags
    d               = 0-Youtvec;
    eq.flag.exist   = ~all(sign(d)==max(sign(d))); 
    eq.flag.equi    = (abs(0-Youtvec(tt))<options.tolY);
    eq.flag.down    = (0<Youtvec(tt));
    eq.flag.up      = (0>Youtvec(tt));
    % 5. Shift bounds
    Ylb             = (eq.flag.up)*Y    + (eq.flag.down)*Ylb;
    Yub             = (eq.flag.up)*Yub  + (eq.flag.down)*Y;
    % 6. Break if equilibrium
    if eq.flag.equi,break,end
    % 7. Plot option
%     if strcmp(options.eqplot,'Y') 
%        figure(888);       
%        subplot(2,2,3);
%        plot(Yinvec(Yinvec~=0),'bo','markersize',6,'markerfacecolor','b');hold on;grid on;
%        plot(Youtvec(Youtvec~=0),'ro','markersize',6,'markerfacecolor','r');
%        xlim([0,max(tt,12)]);
%        legend('Yin','Yout','Location','NorthEast');
%        title('Equilibrium');
%        drawnow;  
%     end
end



end