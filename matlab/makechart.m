function makechart(titlelist,xvalues,legendlist,figlabel,zdata1,zdata2,zdata3)



figure(1)

titlelist = char(strrep(cellstr(titlelist),'_','.'));

ndsets=3;       % default, changed below as applicable
if nargin==5
    zdata2=nan*zdata1;
    zdata3=nan*zdata1;
    ndsets =1;
elseif nargin == 6
    zdata3 =nan*zdata1;
    ndsets=2;
elseif ((nargin>7) | (nargin <=3))
    error ('makechart takes 5 to 7 arguments')
end

nvars = size(titlelist,1);
if nvars==1 
    nrows=1;
    ncols = 1;
elseif nvars==2
    nrows =2;
    ncols = 1;
elseif (nvars == 3 | nvars ==4)
    nrows = 2;
    ncols =2;
elseif (nvars==5 |nvars ==6)
    nrows = 2;
    ncols = 3;
elseif (nvars==7 | nvars==8)
    nrows = 4;
    ncols = 2;
else 
    error('too many variables (makechart)')
end

k = zeros(nvars,1);
for i = 1:nvars
    subplot(nrows,ncols,i)
    h1=plot(xvalues,zdata1(:,i),'k-',xvalues,zdata2(:,i),'r--',xvalues,zdata3(:,i),'b:');
    [x0 x1 y10 y11] = pickaxes(xvalues,zdata1(:,i));
    [x0 x1 y20 y21] = pickaxes(xvalues,zdata2(:,i));
    [x0 x1 y30 y31] = pickaxes(xvalues,zdata3(:,i));
    y0 = min([y10,y20,y30]);
    y1 = max([y11,y21,y31]);
    if y0==y1
        y1=y0+1;
    end
    %y0 = floor(y0*10)/10;
    %y1 = ceil(y1*10)/10;
    if i==1
        y1 = y1 + 0.1;
    end
    
    %x0 = 5+5.*floor((x0-5)./5);
    %x1 = 5+5.*ceil((x1-5)./5);
    
    if (y1-y0 < 0.05)
        y0 = y0-0.001;
        y1 = y1+0.001;
    elseif (y1-y0 <= 0.1)
        y0 = 0.02+0.02*floor((y0-0.02)./0.02);
        y1 = 0.02+0.02.*ceil((y1-0.02)./0.02);
        yt = [y0:0.02:y1];
        set(gca,'YTick',yt);  
    elseif (y1-y0 <= 0.5)
        k(i) = 4;
        y0 = 0.1+0.1*floor((y0-0.1)./0.1);
        y1 = 0.1+0.1.*ceil((y1-0.1)./0.1);
        yt = [y0:0.1:y1];
        set(gca,'YTick',yt);    
    elseif (y1-y0 <= 1)
        k(i) = 1;
        y0 = 0.2+0.2*floor((y0-0.2)./0.2);
        y1 = 0.2+0.2.*ceil((y1-0.2)./0.2);
        yt = [y0:0.2:y1];
        set(gca,'YTick',yt);
    elseif (y1-y0 <= 2.5)
        k(i) = 2;
        y0 = 0.5+0.5*floor((y0-0.5)./0.5);
        y1 = 0.5+0.5.*ceil((y1-0.5)./0.5);
        yt = [y0:0.5:y1];
        set(gca,'YTick',yt);
    elseif (y1-y0 <= 5)
        k(i) = 3;
        y0 = floor(y0);
        y1 = ceil(y1);
        yt = [y0:1:y1];
        set(gca,'YTick',yt);
    else
        y0 = 2+2*floor((y0-2)./2);
        y1 = 2+2.*ceil((y1-2)./2);
        yt = [y0:2:y1];
        set(gca,'YTick',yt);
    end
    
    k    
    %if i==1
    %    y1 = y1 + 0.2;
    %end
    axis([x0 x1 y0 y1])
    set(h1,'linewidth',2);
    
    
    
    xlabel('Quarters');
    %xlabel('Draws');
    if i==1
        if isempty(legendlist) ~= 1
          legend(legendlist,'Location','NorthWest')
          legend('boxoff');
          %ylim([y0 y1*1.1]);
        end  
        text('String',figlabel,'Units','normalized','Position',[1.2 1.24],...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
       
    end
    
    %ylabel('Percent');
    title(strtrim(titlelist(i,:)),'FontSize',11);
    
end

