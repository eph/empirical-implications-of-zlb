%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot_nonlinear_irf.m:  plot impulse responses from nonlinear and
%linear models 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;
close all;

fig_dir = '../';

Neulers = 18;
Nvar = 30;
TT = 20;
compareswitch = 0;  %if compareswitch = 0, compare nonlinear
                    %constrained to linear soluton
                    %if compareswitch = 1, compare nonlinear
                    %constrained to unconstrained, nonlinear
                    %if compareswitch = -1, then no comparison

%shocktype = 1;  %1 is for liquidity shock; 2 is MEI shock
dset = 'mean'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read data from disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(fig_dir,'results/irf/nonlinearirf_', dset, '_',num2str(shocktype),'.txt');
endogirf = load(file);

file = strcat(fig_dir,'results/irf/eulers_', dset, '_',num2str(shocktype),'.txt');
eulersirf = load(file);

if (compareswitch == 0)
    file = strcat(fig_dir,'results/irf/linearirf_', dset, '_',num2str(shocktype),'.txt');
end

if (compareswitch == 1)
    file = strcat(fig_dir,'results/irf/nonlinearirf_unc_', dset, '_',num2str(shocktype),'.txt');
end

if (compareswitch == -1) 
    endogirf2 = NaN*ones(Nvar*TT,1);
else
    endogirf2 = load(file);
end 

endogirf = reshape(endogirf,Nvar,TT);
endogirf2 = reshape(endogirf2,Nvar,TT);
eulersirf = reshape(eulersirf,2*Neulers,TT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


varswitch = 1;
if (varswitch == 0)
    titlelist = char('Nominal Rate','Notional Rate','GDP','Inflation',...
               'Investment','Consumption','Utilization','Price of Investment');
    wishlist = [30,31,7,6,3,2,13,11];
    wishlist = [5,9,7,6,3,2,13,11];
    endogplot = endogirf(wishlist,:)';
    endogplot2 = endogirf2(wishlist,:)';
    
else
    titlelist = char('Nominal Rate','Inflation','GDP','Hours Worked',...
               'Investment','Consumption','Price of Installed Capital', ...
                     'Investment Euler Error');
    wishlist = [5,6,7,12,3,2,11];
    eul_index = [9];
    endogplot = [endogirf(wishlist,:)',eulersirf(eul_index,:)'];
    endogplot2 = [endogirf2(wishlist,:)',eulersirf(Neulers+eul_index,:)']; 
    
    %put inflation and nominal rate at annual rate
    endogplot(:,1:2) = 4*endogplot(:,1:2);
    endogplot2(:,1:2) = 4*endogplot2(:,1:2);


end 

figlabel = [];

if (compareswitch == 0)
    legendlist = char('Constrained Nonlinear','Unconstrained Linear');
end

if (compareswitch == 1)
    legendlist = char('Constrained Nonlinear','Unconstrained Nonlinear');
end

if (compareswitch == -1)
    legendlist = char('Constrained Nonlinear');
end

makechart(titlelist,0:TT-1,legendlist,figlabel,endogplot,endogplot2);

set(figure(1),'PaperPositionMode', 'manual') 
set(figure(1), 'PaperUnits', 'inches');
set(figure(1), 'PaperSize',[8.5 11]);
set(figure(1), 'PaperOrientation','portrait');
set(figure(1), 'PaperPosition', [0.25 0.25 [8.5 11]-0.5]);

%if (shocktype == 1)                     % 
%    print -depsc2 ../figure_drafts/nonlinear_irf_liqshock_1218.eps
%    print -dpdf ../figure_drafts/nonlinear_irf_liqshock_1218.pdf    
%else
%    print -depsc2 ../figure_drafts/nonlinear_irf_meishock_1218.eps
%    print -dpdf ../figure_drafts/nonlinear_irf_meishock_1218.pdf    
%end


%figure;
%makechart(titlelist,0:TT-1,legendlist,figlabel,e,endogplot2);

