%Generates Figure with two subplots
% Top Plot: "Distribution of the Probability of Hitting the ZLB"
% Bottom Plot: "Histogram of Duration of ZLB Spells"
clear all;

%bounds for confidence intervals
lb = 0.025;
ub = 0.975; %.975

fig_switch9 = 0;  %if fig_switch = 0, then don't generate figures
figprint_switch = 0;

fig_switch10 = 0;
fig_switch12 = 0;

ndatasets = 1000; % 1000; testing
nobsdata = 125;  %should be 127
capt = 1000000;
nsamples = floor(capt/nobsdata);
%nsamples = 7874;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data from disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zlbsampledata = zeros(ndatasets*nsamples,1);
duration_spell = [];
ss_count = 0;
for k = 1:ndatasets
    file = strcat('./zlbstat-results/nzlbspells',num2str(k-1,'%04d'),'.txt');
    nspells = load(file);
    
    if (nspells > 0) %case where data exists
        %ss_count = ss_count + 1;
       
       file = strcat('./zlbstat-results/zlbduration',num2str(k-1,'%04d'),'.txt');
       zlbduration = load(file);
       %duration_spell = [zlbduration; zlbduration];
       duration_spell = [duration_spell; zlbduration];
       file = strcat('./zlbstat-results/zlbfrequency',num2str(k-1,'%04d'),'.txt');
       zlbfrequency = load(file); 
       %zlbsampledata((ss_count-1)*nsamples+1:nsamples*ss_count) = zlbfrequency;  
       zlbsampledata((ss_count+1):(ss_count+nsamples)) = zlbfrequency;  
       ss_count = ss_count + nsamples; 
    end     
end

fprintf('ss_count = %4.0f\n',ss_count);

zlbsampledata = zlbsampledata/100;
zlbsampledata = zlbsampledata(1:ss_count); %-1); % error here bc of limited size - added -1

%% Right Tail of Sample Histogram

ind_zerozlb = find(zlbsampledata <= 1e-08);
zlbzerofreq = length(ind_zerozlb)/length(zlbsampledata);
zlbfreqmean = mean(zlbsampledata);
zlbfreqmax = max(zlbsampledata);
zlbfreqss_s = sort(zlbsampledata);
 
percentile = [0.95; 0.975; 0.995];
for i = 1:3
    ordinal_rank = round(percentile(i)*length(zlbsampledata));
    zlbsamplefreq_upper_ci(i) = zlbfreqss_s(ordinal_rank);
end


Ndur = length(duration_spell);
duration_spell_s = sort(duration_spell);
durationspellmax = max(duration_spell);

ind_short = find(duration_spell <= 4);
ind_long = find(duration_spell >= 12);
freq_short = length(ind_short)/Ndur;
freq_long = length(ind_long)/Ndur;
freq_med = 1-freq_short-freq_long;

percentile = [0.95; 0.975; 0.995];
for i = 1:3
    ordinal_rank_dur = round(percentile(i)*length(duration_spell_s)+1);
    zlbdur_upper_ci(i) = duration_spell_s(ordinal_rank_dur);
end





  figure1a = figure(1);
  fw = 'normal';
    
        %handaxes1 = axes('position', [0.1306 0.1077 0.7787 0.226]);
        subplot(2,1,1);
        
        edges4a = [5:2:50];
        Ndata = length(zlbsampledata);
        [Nzlb4a] = histc(100*zlbsampledata,edges4a);
        freq4a = 100*Nzlb4a./Ndata;
        hb=bar(edges4a,freq4a,'style','histc');
        set(hb,'FaceColor','b')
        colormap([.15 .15 .15])
                 title('Right Tail of Sample Histogram of Prob($0\leq R_t \leq \overline{R}$)',...
             'interpreter','latex','FontName',...
             'times','FontWeight',fw,'Fontsize',14);
        %title('Right Tail of Histogram of Duration of ZLB Spells','FontName',...
        %      'times','FontWeight',fw,'Fontsize',14);
        ylabel('Percent of Observations','FontName',...
            'times','FontWeight',fw,'Fontsize',13);
        xlabel('Percent','FontName',...
            'times','FontWeight',fw,'Fontsize',13);
        xlim([2 30]);
        
        hold on;
        %plot([zlbsamplefreq_upper_ci(1)*100 zlbsamplefreq_upper_ci(1)*100],[0 2],'k:');
        %plot([zlbsamplefreq_upper_ci(2)*100 zlbsamplefreq_upper_ci(2)*100],[0 4],'k:');
        hold off;
        %text(zlbsamplefreq_upper_ci(1)*100 - 3,2.15,'97.5th Percentile','FontName','Times New Roman','FontWeight','bold','FontSize',12);
        %text(zlbsamplefreq_upper_ci(2)*100 - 7,4.15,'95th Percentile','FontName','Times New Roman','FontWeight','bold','FontSize',12);
        %text(15,9,'Average = ','FontName','Times New Roman','FontWeight','bold','FontSize',12);
        %text(20,9,sprintf('%9.1f%%',100*zlbfreqmean),'FontName','Times New Roman','FontWeight','bold','FontSize',12);
        xlim([5 50]);
        %axis([2 30 0 12]);
        %set(handaxes1, 'box', 'off');
        box off;
        
        
      
        
        
        % Duration Figure
        subplot(2,1,2);
        
        edges1 = [5:2:40];
        [Ndurbin1] = histc(duration_spell,edges1);
        freq1 = 100*Ndurbin1./Ndur;
        hb=bar(edges1,freq1,'histc');
        set(hb,'FaceColor','b')
        title('Right Tail of Histogram of Duration of ZLB Spells','FontName',...
              'times','FontWeight',fw,'Fontsize',14);
        ylabel('Percent of Spells','FontName',...
               'times','FontWeight',fw,'Fontsize',13);
        xlabel('Quarters','FontName',...
               'times','FontWeight',fw,'Fontsize',13);
        %axis([5 30 0 8]);
        xlim([5 40]);
        %set(handaxes1, 'box', 'off');
        box off;
        hold on;
        %plot([zlbdur_upper_ci(1) zlbdur_upper_ci(1)],[0 2],'k:');
        %plot([zlbdur_upper_ci(2) zlbdur_upper_ci(2)],[0 1],'k:');
        hold off;
        %text(zlbdur_upper_ci(1) - .5,2.3,'97.5th Percentile','FontName','Times New Roman','FontWeight','bold','FontSize',12);
        %text(zlbdur_upper_ci(2) - .5,1.3,'99.5th Percentile','FontName','Times New Roman','FontWeight','bold','FontSize',12);

        
        % upper inner plot
        
           handaxes2 = axes('position', [0.5502 0.7492 0.3268 0.1351]);   
        edges2 = [20:5:70];
        [Nzlbbin2] = histc(100*zlbsampledata,edges2);
        freq2 = 100*Nzlbbin2./Ndata;
        hb=bar(edges2,freq2,'histc');   
        set(hb,'FaceColor','b')
        title('Far-Right Tail of the Histogram','FontName',...
              'times','interpreter','latex','FontWeight',fw,'Fontsize',10);
        ylabel('Percent of Observations','FontName',...
               'times','FontWeight',fw,'Fontsize',10);
        xlabel('Percent','FontName',...
               'times','FontWeight',fw,'Fontsize',10);
        xlim([20 70]);
        set(handaxes2,'XTick',[20:10:70]);
        set(handaxes2, 'box', 'off');
        
        % lower inner plot 
         handaxes2 = axes('position', [0.5502 0.2692 0.3268 0.1351]);   
        edges3 = [20:5:70];
        [Ndurbin3] = histc(duration_spell,edges3);
        freq3 = 100*Ndurbin3./Ndur;
        hb=bar(edges3,freq3,'histc');   
        set(hb,'FaceColor','b')
        title('Far-Right Tail of the Histogram','FontName',...
              'times','interpreter','latex','FontWeight',fw,'Fontsize',10);
        ylabel('Percent of Spells','FontName',...
               'times','FontWeight',fw,'Fontsize',10);
        xlabel('Quarters','FontName',...
               'times','FontWeight',fw,'Fontsize',10);
        xlim([20 70]);
        %ylim([0 0.15]);
        %set(handaxes2,'YTick',[0 0.05 0.1 0.15]);
        set(handaxes2,'XTick',[20:10:70]);
        set(handaxes2, 'box', 'off');
        
        
        set(figure(1),'PaperPositionMode', 'manual') 
        set(figure(1), 'PaperUnits', 'inches');
        set(figure(1), 'PaperSize', [6 8]);
        %set(figure(1), 'PaperPosition', [0.05 0.05 5.9 7.9]); 
        set(figure(1), 'PaperPosition', [0 -0.4 6 8.5]); 
        %print -dpsc2 'figure_drafts/zlbfreqhist.ps';
        %print -dpdf 'figure_drafts/zlbfreqhist.pdf';
        %print -dpdf '/msu/res2/Shared_Projects/GLS/G_LS_S/latex_aer_revision/zlb_freq_duration.pdf';entile','FontName','Times New Roman','FontWeight','bold','FontSize',12);
        %text(zlbdur_upper_ci(2) - .5,1.3,'99.5th Percentile','FontName','Times New Roman','FontWeight','bold','FontSize',12);


        
        
        % upper inner plot
%         
%            handaxes2 = axes('position', [0.5502 0.7492 0.3268 0.1351]);   
%         edges2 = [30:5:70];
%         [Nzlbbin2] = histc(100*zlbsampledata,edges2);
%         freq2 = 100*Nzlbbin2./Ndata;
%         hb=bar(edges2,freq2,'histc');   
%         set(hb,'FaceColor','b')
%         title('Far-Right Tail of the Histogram','FontName',...
%               'times','interpreter','latex','FontWeight',fw,'Fontsize',10);
%         ylabel('Percent of Observations','FontName',...
%                'times','FontWeight',fw,'Fontsize',10);
%         xlabel('Percent','FontName',...
%                'times','FontWeight',fw,'Fontsize',10);
%         xlim([30 70]);
%         set(handaxes2, 'box', 'off');
%         
%         % lower inner plot 
%          handaxes2 = axes('position', [0.5502 0.2692 0.3268 0.1351]);   
%         edges3 = [20:5:70];
%         [Ndurbin3] = histc(duration_spell,edges3);
%         freq3 = 100*Ndurbin3./Ndur;
%         hb=bar(edges3,freq3,'histc');   
%         set(hb,'FaceColor','b')
%         title('Far-Right Tail of the Histogram','FontName',...
%               'times','interpreter','latex','FontWeight',fw,'Fontsize',10);
%         ylabel('Percent of Spells','FontName',...
%                'times','FontWeight',fw,'Fontsize',10);
%         xlabel('Quarters','FontName',...
%                'times','FontWeight',fw,'Fontsize',10);
%         xlim([20 70]);
%         ylim([0 0.15]);
%         set(handaxes2,'YTick',[0 0.05 0.1 0.15]);
%         set(handaxes2, 'box', 'off');
%         
        
%         set(figure(1),'PaperPositionMode', 'manual') 
%         set(figure(1), 'PaperUnits', 'inches');
%         set(figure(1), 'PaperSize', [6 8]);
%         set(figure(1), 'PaperPosition', [0.05 0.05 5.9 7.9]); 
        %print -dpsc2 'figure_drafts/zlbfreqhist.ps';
        %print -dpdf 'figure_drafts/zlbfreqhist.pdf';
        
        print -dpdf '/msu/res2/Shared_Projects/GLS/G_LS_S/latex_aer_revision/figures_aer2/zlb_freq_duration_rbar25_1000draws.pdf'; %%% print command used
        %saveas(figure1a,'/msu/res2/Shared_Projects/GLS/G_LS_S/latex_aer_revision/figures_aer2/zlb_freq_duration_rbar10_1000draws.pdf','pdf');
        
        % 508
%         upper_main = [edges4a' freq4a];
%         upper_inner = [edges2' freq2];
%         lower_main = [edges1' freq1];
%         lower_inner = [edges3' freq3];
%         xlswrite('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/zlbfreqhist_upper_main.xls',upper_main);
%         xlswrite('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/zlbfreqhist_upper_inner.xls',upper_inner);
%         xlswrite('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/zlbfreqhist_lower_main.xls',lower_main);
%         xlswrite('/msu/res5/Shared_Projects/GHLSS_AER/FEDS/508/zlbfreqhist_lower_inner.xls',lower_inner);
%        