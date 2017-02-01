%Generates Figure "Distribution of the Probability of Hitting the ZLB"
clear all;

capt = 1000;
%ndraws = 1000;
ndraws = 1000;
nmom = 4;
modelswitch = 1; %0 for nonlinear,constrained; 1 for nonlinear,
                 %unconstrained; 2 for linear


filterswitch = 0;  %0 for hpfilter; 1 for bkfilter; 2 for cffilter

Nvar = 34;
TT = capt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data from disk and calculate selected moments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndata = 125;
%obs(1)=dgdp,obs(2)=dc,obs(3)=dinv,obs(4)=dh,obs(5)=dp,obs(6)=ff'
% nobs = 7;
% usdata = load('glss_data.txt');
% 
% %compute log-level of output, consumption and investment
% nwish = nmom;
% usexp = zeros(ndata+1,nwish);
% for i = 1:nwish
%     for tt = 2:ndata+1
%         usexp(tt,i) = usexp(tt-1,i) + usdata(tt-1,i);
%     end    
% end
% 
% usexp_level = usexp(2:ndata+1,:);
% stdvec_usdata = zeros(nwish,1);
% for i = 1:nwish
%     if (filterswitch == 0)
%         ydev = hpfast(usexp_level(:,i),1600);
%     elseif (filterswitch == 1)
%         ydev = bkfilter(usexp_level(:,i),6,32);
%     else
%         ydev = cffilter(usexp_level(:,i),6,32);
%     end
%     stdvec_usdata(i) = std(ydev);
% end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data from disk and calculate selected moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stdvecall = zeros(nmom,ndraws);
countdsets = 0;
for k = 1:ndraws
    
    %parameter values from disk
    file = strcat('./final-final/by-draw/parasim',num2str(k-1,'%04d'),'.txt');
    paramdata = load(file);
    
    if (modelswitch == 0) 
        file = strcat('./zlbstat-results/modeldata_file',num2str(k-1,'%04d'),'.txt');
    elseif (modelswitch == 1)
        file = strcat('./zlbstat-results/modeldata_unc_file',num2str(k-1,'%04d'),'.txt');
    else
        file = strcat('./zlbstat-results/modeldata_linear_file',num2str(k-1,'%04d'),'.txt');
    end
    
    
    modeldata = load(file);
    if (length(modeldata) > 1) 
        
        modeldata = reshape(modeldata,Nvar,TT)';
        gz = paramdata(3,1);

        techshk = modeldata(:,25);
        tech_level = zeros(capt+1,1);
        for i = 2:capt+1
            tech_level(i) = log(gz) + tech_level(i-1) + techshk(i-1);
        end
        tech_level = tech_level(2:capt+1);

        wishlist = [7,2,3];  %gdp,cons,inv
        nwish = nmom-1; 
        stdvec = zeros(nmom,1);
        for i = 1:nwish
            data = log(modeldata(:,wishlist(i))) + tech_level;
            data = data(capt-ndata+1:capt);
            if (filterswitch == 0)
                ydev = hpfast(data,1600);
            elseif (filterswitch == 1)
                ydev = bkfilter(data,6,32);
            else
                ydev = cffilter(data,6,32);
            end
            stdvec(i) = std(ydev);
        end

        labor = log(modeldata(capt-ndata+1:capt,12));
        if (filterswitch == 0)
            ydev = hpfast(labor,1600);
        elseif (filterswitch == 1)
            labor = labor(1001:11000);
            ydev = bkfilter(labor,6,32);
        else
            labor = labor(1001:11000);
            ydev = cffilter(labor,6,32);
        end
        stdvec(4) = std(ydev);
        countdsets = countdsets + 1;
        stdvecall(:,countdsets) = stdvec;
    end 
end 

stdvecall2 = stdvecall(:,1:countdsets);
stdvecavg = mean(stdvecall2,2);
varavg = mean(stdvecall2.^2,2);
stdvecavg2 = sqrt(varavg);









