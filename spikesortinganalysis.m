function [Ir, Xc, tvec, spikes, spikes_sorted, spikes_sorted2, VS_max, VS_max2, comments] = spikesortinganalysis(datapath,manualyn,plotyn,saveyn,savefigyn)
%
% ----------------------
% Spike Sorting Analysis
% ----------------------
%
% For extracellular recordings of neurons, this program employs both
% principal components analysis (PCA) and wavelet theory to extract the
% features of all spikes in a time series. Clustering is performed with a 
% k-means algorithm, and the optimal number of clusters with the gap method.
% It then determines the degree of entrainment with a stimulus by calculating 
% vector strength for each spike cluster. Other metrics are additionally
% calculated for all sorted spikes.
%
% [Ir, Xc, tvec, spikes, spikes_sorted, spikes_sorted2, VS_max, VS_max2] =
%       spikesortinganalysis(datapath,manualyn,plotyn,saveyn,savefigyn)
%
%
% INPUTS:
%   - datapath: The folder in which an "Extracted Data.mat" file was
%   obtained using the function importVICdata()
%   - manualyn: If set to 1, will ask the user to manually define an
%   optimal number of clusters. (default = 0)
%   - plotyn: If set to 1, will plot various metrics from both sorting
%   algorithms. (default = 0)
%   - saveyn: If set to 1, will save the data in a new file. (default = 0)
%   - savefigyn: If set to 1, will save the figures in a new file. (default = 0)
%
% OUTPUTS:
%   - spikes: All of the unsorted spikes, separated by time series.
%   - spikes_sorted(2): sorted spikes from PCA and wavelet theory
%   - VS_max(2): maximum vector strength across clusters for PCA and
%   wavelet theory
%
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Loading and Initialization       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
close all;clc


if datapath(end) == '/'
    disp(['Loading ' datapath 'Extracted Data.mat...']);
    load([datapath 'Extracted Data.mat']);
else
    disp(['Loading ' datapath '/Extracted Data.mat...']);
    load([datapath '/Extracted Data.mat']);
end

if exist('manualyn') == 0
    manualyn = 0;
end
if exist('plotyn') == 0
    plotyn = 0;
end
if exist('saveyn') == 0
    saveyn = 0;
end
if exist('savefigyn') == 0
    savefigyn = 0;
end



disp('Sorting and analyzing spikes...');

% Spike analysis
setfiguredefaults(10);
tel = 0;
tel1 = 0;
tel2 = 0;

% If using ABF data
nonraw = 1:length(Ir);

% initialization
for m = 1:length(nonraw)
    tic;
    
        tel = tel1 + tel2 + tel;
        disp(['Iteration ' num2str(m) ' out of ' num2str(length(nonraw)) '; Total time elapsed: ' num2str(tel) ' seconds.']);
    
        % Define command signal and neural response
        xc = Xc{nonraw(m)};
        x = Ir{nonraw(m)};
        
        % Calculate entropy (normalized by log(N) of stimulus)
        entropy_stim(m) = entropy1(xc,2^8)/log10(2^8);
        
        % Sample rate (Hz)
        sr = Fs;
        
        % Spike length for extraction (s)
        sl = 0.004;
        sl = sl*sr;
        
        % Extract spikes, peaks in command signal, set maximum frequency
        % limit for optimization purposes
        thr = 6*median(abs(x)./0.6745);
        [pk{m}, tr{m}] = PTDetect(x, thr);
        [pkstim{m}, trstim{m}] = PTDetect(xc, std(xc-mean(xc)));
        numpkstim = length(pkstim{m});
        sl_stim = length(xc)/numpkstim;
        stim_freq(m) = numpkstim/length(xc)*sr;
        max_freq = 1e3;
        
    
if ~isempty(tr{m}) && length(tr{m}) < 5e4
    for j = 1:length(tr{m})  
        
            % Extract spikes
            try
                spikes{m}(:,j) = x(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2));
                spikes_stim{m}(:,j) = xc(tr{m}(j) - round(sl*12.5):tr{m}(j) + round(sl*12.5)); 
            catch
                spikes{m}(:,j) = zeros(1,length(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2)));
                spikes_stim{m}(:,j) = zeros(1,length(tr{m}(j) - round(sl*12.5):tr{m}(j) + round(sl*12.5)));                                 
            end      
            
            % Extract spike amplitudes, power, and areas
            try
                pkpkampl{m}(j) = max(spikes{m}(:,j)) - min(spikes{m}(:,j));
                spike_power{m}(j) = trapz(tvec{m}(1:length(spikes{m}(:,j))),spikes{m}(:,j).^2);
                spike_area{m}(j) = trapz(tvec{m}(1:length(spikes{m}(:,j))),spikes{m}(:,j));
            catch
                disp('Error. Line 132.');
            end
            
            if j > 1
                try
                    IEI{m}(j-1) = (tr{m}(j) - tr{m}(j-1))/sr;
                end
            end
            
                  
    end
    try
        IEI_mean(m) = mean(IEI{m});
        IEI_std(m) = std(IEI{m});
        IEI_CV(m) = IEI_std(m)/IEI_mean(m);
    end
        try
            % Mean amplitude (all spikes)
            pkpkampl_mean(m) = mean(pkpkampl{m});
            pkpkampl_sem(m) = std(pkpkampl{m})/sqrt(length(tr{m}));
        end
    
        try
            % Mean power (all spikes)
            spike_power_mean(m) = mean(spike_power{m});
            spike_power_sem(m) = std(spike_power{m})/sqrt(length(tr{m}));
        end

        try
            % Mean area (all spikes)
            spike_area_mean(m) = mean(spike_area{m});
            spike_area_sem(m) = std(spike_area{m})/sqrt(length(tr{m}));
        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting by Principal components analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = 1;
    clear nc dim1
    
    try
        nc = numclust(m);
    end
    try 
        dim1 = dims(m);
    end
    m2=1;
    if exist('dim1') == 0
        dims(m) = 2;
    end
    nn=0;
     while r == 1 && exist('nc') == 0 && nn < 4
         %}
        
        [c{m}, score{m}, latent{m}, tsquared{m}, explained{m}, mu{m}]=pca(spikes{m});
        
        close all;figure;
        try
            subplot(2,1,1);plot(x);hold on;plot(xc.*(range(x)/range(xc).*0.5),'r');
            subplot(2,1,2);scatter(c{m}(:,1),c{m}(:,2));
            xlabel(['PC1 (' num2str(explained{m}(1)) '%)']);ylabel(['PC2 (' num2str(explained{m}(2)) '%)']);
        end
        
        clear clust1
        
        % Explained variance should be at least 80% across dims(m)
        % principal components.
        explainedvariance(m) = 0;
        nn=0;
        while explainedvariance(m) < 70 && nn < 4
            
        if manualyn == 1
            if m2 > 1
                numclust(m) = input('Number of clusters:   ');
                dims(m) = input('Dimensions (1,2,3,...):   ');
            else
                evalclust{m} = evalclusters(c{m}(:,1:dims(m)),'gmdistribution','gap','KList',1:4);
                numclust(m) = evalclust{m}.OptimalK;
            end
        else
            try
              
                evalclust{m} = evalclusters(c{m}(:,1:dims(m)),'gmdistribution','gap','KList',1:4);

                numclust(m) = evalclust{m}.OptimalK;
            end
        end
        try
            GMmodel{m} = fitgmdist(c{m}(:,1:dims(m)),numclust(m),'Replicates',5);
            T{m} = cluster(GMmodel{m},c{m}(:,1:dims(m)));
            explainedvariance(m) = sum(explained{m}(1:dims(m)));
        catch 
            T{m} = [];
            explainedvariance(m) = 0;
        end
        
        if explainedvariance(m) < 70
            dims(m) = dims(m) + 1;
        end
        nn=nn+100;
        end
        
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r')
        try
        if numclust(m) >= 2
            
        for k = 1:numclust(m)
            subplot(2,1,2);scatter(c{m}(T{m}==k,1),c{m}(T{m}==k,2));
            title([num2str(numclust(m)) ' clusters'])
            xlabel(['PC1 (' num2str(explained{m}(1)) '%)']);ylabel(['PC2 (' num2str(explained{m}(2)) '%)']);
            hold all;
        end
        
        else
            k=1;subplot(2,1,2);scatter(c{m}(:,1),c{m}(:,2));
        end
        end
        if manualyn == 1
            r = input('Repeat? (1=yes):   ');
        else
            r = 0;
        end
        m2=m2+1;
        nn=nn+1;
    end
    
    
    try
        [c{m}, score{m}, latent{m}, tsquared{m}, explained{m}, mu{m}]=pca(spikes{m});
%         T{m} = kmeans(c{m}(:,1:dims(m)),numclust(m),'Replicates',5);
        GMmodel{m} = fitgmdist(c{m}(:,1:dims(m)),numclust(m),'Replicates',5);
        T{m} = cluster(GMmodel{m},c{m}(:,1:dims(m)));
    catch
        disp('Unable to find number of clusters. Line 246.');
    end
%}
    
    % Sort the spikes and calculate statistics for each subset
    try
    if numclust(m) >=2
    for k = 1:numclust(m)
        spikes_sorted{m}{k} = spikes{m}(:,T{m}==k);
        spikes_stim_sorted{m}{k} = spikes_stim{m}(:,T{m}==k);
        entropy_spikes_stim_sorted(m,k) = entropy1(mean(spikes_stim_sorted{m}{k},2),2^8)/log10(2^8);
        entropy_diff(m,k) = entropy_spikes_stim_sorted(m,k) - entropy_stim(m);
        tr_sorted{m}{k} = tr{m}(T{m}==k);
%         pk_sorted{m}{k} = pk{m}(T{m}==k);
        for j = 2:length(tr_sorted{m}{k})
            IEI_sorted{m}{k}(j-1)=(tr_sorted{m}{k}(j)-tr_sorted{m}{k}(j-1))/sr;
        end
        try
            IEI_sorted_mean(m,k) = mean(IEI_sorted{m}{k});
            IEI_sorted_std(m,k) = std(IEI_sorted{m}{k});
            IEI_sorted_CV(m,k) = IEI_sorted_std(m,k)/IEI_sorted_mean(m,k);
        catch
            disp('Error. Line 272.');
        end
    end
    else
        spikes_sorted{m}{1} = spikes{m};
        spikes_stim_sorted{m}{k} = spikes_stim{m};
        tr_sorted{m}{1} = tr{m};
        for j = 2:length(tr_sorted{m}{1})
            IEI_sorted{m}{1}(j-1)=(tr_sorted{m}{1}(j)-tr_sorted{m}{1}(j-1))/sr;
        end
%         pk_sorted{m}{1} = pk{m};
        try
            IEI_sorted_mean(m,k) = mean(IEI_sorted{m}{k});
            IEI_sorted_std(m,k) = std(IEI_sorted{m}{k});
            IEI_sorted_CV(m,k) = IEI_sorted_std(m,k)/IEI_sorted_mean(m,k);
        catch
            disp('Error. Line 288.');
        end
    end
    end
    
   
    try
    for k = 1:numclust(m)
        try
        spikes_sorted_pow_mean{m}(k) = mean(spike_power{m}(T{m}==k));
        spikes_sorted_pow_sem{m}(k) = std(spike_power{m}(T{m}==k))/sqrt(sum(T{m}==k));
        spikes_sorted_area_mean{m}(k) = mean(spike_area{m}(T{m}==k));
        spikes_sorted_area_sem{m}(k) = std(spike_area{m}(T{m}==k))/sqrt(sum(T{m}==k));
        spikes_sorted_pkpkampl_mean{m}(k) = mean(pkpkampl{m}(T{m}==k));
        spikes_sorted_pkpkampl_sem{m}(k) = std(pkpkampl{m}(T{m}==k))/sqrt(sum(T{m}==k));
        end
    end
    end
    
    tel1 = toc;
    try
        disp(['PCA: Number of Clusters = ' num2str(numclust(m)) '; Time elapsed: ' num2str(tel1) ' seconds.']);  
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort spikes using the WAVELET transform %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

for k = 1:size(spikes{m},2)
    % Discrete time wavelet transform
    [cA{m}(:,k),cD{m}(:,k)] = dwt(spikes{m}(:,k),'haar');
end
for l = 1:size(cA{m},1)
    % Find the coefficients that show greatest spread using Lilliefors test
    % A large statistic corresponds to a divergence from a normal
    % distribution
    warning off
    try
    [~,~,kstat{m}(l)] = lillietest(cA{m}(l,:));
    end
end

% Sort by K-statistic
try
sortedK = sort(kstat{m});
coeff1(m) = find(kstat{m} == sortedK(end));
end
try
coeff2(m) = find(kstat{m} == sortedK(end-1));
end
try
coeff3(m) = find(kstat{m} == sortedK(end-2));
end

% Cluster data
clear nc
r = 1;
    try
        nc = numclust2(m);
    end
    try 
        dim1 = dims(m);
    end
m2=1;
    if exist('dim1') == 0
        dims(m) = 2;
    end
     while r == 1 && exist('nc') == 0 
         %}
        
       try
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r');
        subplot(2,1,2);scatter(cA{m}(coeff1(m),:),cA{m}(coeff2(m),:));
        xlabel('C1');ylabel('C2');
       end
        
                
        clear clust1 clust2
        
        if manualyn == 1
            if m2 > 1
                numclust2(m) = input('Number of clusters:   ');
            else
%                 for j = 2:5
%                     clust2(:,j-1) = kmeans([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:)]',j,'Replicates',5);
%                 end
                evalclust2{m} = evalclusters([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]','gmdistribution','gap','KList',1:4);
                numclust2(m) = evalclust2{m}.OptimalK;
            end
        else
%             for j = 2:5
%                 clust2(:,j-1) = kmeans([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:)]',j,'Replicates',5);
%             end
try
            evalclust2{m} = evalclusters([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]','gmdistribution','gap','KList',1:4);
            numclust2(m) = evalclust2{m}.OptimalK;
end
        end
        
        try
%           T2{m} = kmeans([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:)]',numclust2(m),'Replicates',5);
            GMmodel2{m} = fitgmdist([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]',numclust2(m),'Replicates',5);
            T2{m} = cluster(GMmodel2{m},[cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]');
        catch
            disp('Error. Line 370.');
            T2{m} = [];
        end
        
         
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r')
        try
        if numclust2(m) >= 2
        colors = {'r','g','b','y','k','m','c'};    
        for k = 1:numclust2(m)
            subplot(2,1,2);scatter(cA{m}(coeff1(m),T2{m}==k),cA{m}(coeff2(m),T2{m}==k),colors{k});
            title([num2str(numclust2(m)) ' clusters'])
            xlabel('C1');ylabel('C2');
            hold all;
        end
        
        else
            
            k=1;subplot(2,1,2);scatter(cA{m}(coeff1(m),:),cA{m}(coeff2(m),:));
            
        end
        end
        if manualyn == 1
            r = input('Repeat? (1=yes):   ');
        else 
            r=0;
        end
        m2=m2+1;
    
    end
    
    
    try
%         T2{m} = kmeans([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:)]',numclust2(m),'Replicates',5);
        GMmodel2{m} = fitgmdist([cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]',numclust2(m),'Replicates',5);
        T2{m} = cluster(GMmodel2{m},[cA{m}(coeff1(m),:);cA{m}(coeff2(m),:);cA{m}(coeff3(m),:)]');
    catch
        disp('Unable to find number of clusters. Line 407.');
    end
    
% Sort the spikes by clustering from wavelet coefficients
try
    if numclust2(m) >=2
    for k = 1:numclust2(m)
        try
        spikes_sorted2{m}{k} = spikes{m}(:,T2{m}==k);
        end
        try
        spikes_stim_sorted2{m}{k} = spikes_stim{m}(:,T2{m}==k);
        end
        try
        entropy_spikes_stim_sorted2(m,k) = entropy1(mean(spikes_stim_sorted2{m}{k},2),2^8)/log10(2^8);
        entropy_diff2(m,k) = entropy_spikes_stim_sorted2(m,k) - entropy_stim(m);
        end
        try
        tr_sorted2{m}{k} = tr{m}(T2{m}==k);
        end
%       pk_sorted{m}{k} = pk{m}(T{m}==k);
            for j = 2:length(tr_sorted2{m}{k})
                IEI_sorted2{m}{k}(j-1)=(tr_sorted2{m}{k}(j)-tr_sorted2{m}{k}(j-1))/sr;
            end
            try
        IEI_sorted2_mean(m,k) = mean(IEI_sorted2{m}{k});
        IEI_sorted2_std(m,k) = std(IEI_sorted2{m}{k});
        IEI_sorted2_CV(m,k) = IEI_sorted2_std(m,k)/IEI_sorted2_mean(m,k);
            catch
                disp('Error. Line 427.');
            end
    end
    else
        spikes_sorted2{m}{1} = spikes{m};
        spikes_stim_sorted2{m}{k} = spikes_stim{m};
        tr_sorted2{m}{1} = tr{m};
%         pk_sorted{m}{1} = pk{m};
try
        for j = 2:length(tr_sorted2{m}{1})
            IEI_sorted2{m}{1}(j-1)=(tr_sorted2{m}{1}(j)-tr_sorted2{m}{1}(j-1))/sr;
        end
        IEI_sorted2_mean(m,k) = mean(IEI_sorted2{m}{k});
        IEI_sorted2_std(m,k) = std(IEI_sorted2{m}{k});
        IEI_sorted2_CV(m,k) = IEI_sorted2_std(m,k)/IEI_sorted2_mean(m,k);
catch
    disp('Error. Line 443.');
end
    end
end
    try
    for k = 1:numclust2(m)
        try
            spikes_sorted_pow_mean2{m}(k) = mean(spike_power{m}(T2{m}==k));
            spikes_sorted_pow_sem2{m}(k) = std(spike_power{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
            spikes_sorted_area_mean2{m}(k) = mean(spike_area{m}(T2{m}==k));
            spikes_sorted_area_sem2{m}(k) = std(spike_area{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
            spikes_sorted_pkpkampl_mean2{m}(k) = mean(pkpkampl{m}(T2{m}==k));
            spikes_sorted_pkpkampl_sem2{m}(k) = std(pkpkampl{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
        catch
            disp('Error. Line 464.');
        end
    end
    end
   
    tel2 = toc;
    try
    disp(['Wavelet: Number of Clusters = ' num2str(numclust2(m)) '; Time elapsed: ' num2str(tel2) ' seconds.']);
    end
    
    % Extract phases

    % Hilbert transform of command signal
    xc_hilb = hilbert(xc);
    xc_phase = atan2(imag(xc_hilb),real(xc_hilb));
    try
    xc_phase_spikes{m} = xc_phase(tr{m});
    
    % Mean phase and vector strength
    meanphase(m) = mean(xc_phase_spikes{m});
    semphase(m) = std(xc_phase_spikes{m})/sqrt(length(xc_phase_spikes{m}));
    VS(m) = 1/length(tr{m})*abs(sum(exp(sqrt(-1).*xc_phase_spikes{m})));
    
    % Velocity of command signal
    xc_diff = smooth(diff(smooth(Xc{nonraw(m)},10)),10);
    xc_diff_spikes = xc_diff(tr{m});
    xc_diff_spikes_mean(m) = mean(xc_diff_spikes);
    xc_diff_spikes_sem(m) = std(xc_diff_spikes)/sqrt(length(tr{m}));
    catch
        disp('Error. Line 493.')
    end
    
    clear pks_present
    try
    % Bernoulli trials: Determine whether a spike occurs in each cycle of
    % the stimulus.
    if stim_freq(m) < max_freq
    for j = 1:numpkstim
        pks_present(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}))>0;
    end
    end
    catch
        disp('Error. Line 505.');
    end
    
    try
    % Probability of a spike occurring in each cycle of the stimulus
    scaling_fac_VS = sum(pks_present)/numpkstim;
    
    % Scale vector strength by aforementioned probability
    if scaling_fac_VS > 1
        warning(['Scaling factor for vector strength exceeds 1 at index ' num2str(m) '.']);
    elseif scaling_fac_VS < 0
        warning(['Scaling factor for vector strength below 0 at index ' num2str(m) '.']);
    else
        VSsc(m) = VS(m) * scaling_fac_VS;
    end
    catch
        disp('Error. Line 522.');
    end

    
    % Repeat for sorted spikes (PCA)
    try
    for k = 1:numclust(m)
        try
            xc_phase_spikes_sorted{m}{k} = xc_phase_spikes{m}(T{m}==k);
            meanphase_sorted{m}(k) = mean(xc_phase_spikes_sorted{m}{k});
            semphase_sorted{m}(k) = std(xc_phase_spikes_sorted{m}{k})/sqrt(sum(T{m}==k));
            VS_sorted(m,k) = 1/length(tr{m}(T{m}==k))*abs(sum(exp(sqrt(-1).*xc_phase_spikes_sorted{m}{k})));
        catch
            disp('Unable to calculate VS. Line 489.');
        end
        try
        if stim_freq(m) < max_freq
        clear pks_present
        for j = 1:numpkstim
            pks_present(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}(T{m}==k)))>0;
        end
        scaling_fac_VS_sorted(m,k) = sum(pks_present)/numpkstim;
        end
        catch
            disp('Error. Line 500.');
        end
        try
        if scaling_fac_VS_sorted(m,k) > 1
            warning(['Scaling factor for vector strength exceeds 1 at index ' num2str(m) '.']);
        elseif scaling_fac_VS_sorted(m,k) < 0
            warning(['Scaling factor for vector strength below 0 at index ' num2str(m) '.']);
        else
            VSsc_sorted(m,k) = VS_sorted(m,k) * scaling_fac_VS_sorted(m,k);
        end
        catch
            disp('Error. Line 511.');
        end
        
    end
    end
    
    % Repeat for sorted spikes (WAVELET)
    try
    for k = 1:numclust2(m)
        try
            xc_phase_spikes_sorted2{m}{k} = xc_phase_spikes{m}(T2{m}==k);
            meanphase_sorted2{m}(k) = mean(xc_phase_spikes_sorted2{m}{k});
            semphase_sorted2{m}(k) = std(xc_phase_spikes_sorted2{m}{k})/sqrt(sum(T2{m}==k));
            VS_sorted2(m,k) = 1/length(tr{m}(T2{m}==k))*abs(sum(exp(sqrt(-1).*xc_phase_spikes_sorted2{m}{k})));
        catch
            disp('Error. Line 524.');
        end
        try
        if stim_freq(m) < max_freq
        clear pks_present
        for j = 1:numpkstim
            pks_present2(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}(T2{m}==k)))>0;
        end
        end
        catch
            disp('Error. Line 534.');
        end
        
        try
            scaling_fac_VS_sorted2(m,k) = sum(pks_present2(j))/numpkstim;
        catch
            disp('Error. Line 540.');
        end
        
        try
            if scaling_fac_VS_sorted2(m,k) > 1
                warning(['Scaling factor for vector strength exceeds 1 at index ' num2str(m) '.']);
            elseif scaling_fac_VS_sorted(m,k) < 0
                warning(['Scaling factor for vector strength below 0 at index ' num2str(m) '.']);
            else
                VSsc_sorted2(m,k) = VS_sorted2(m,k) * scaling_fac_VS_sorted2(m,k);
            end
        catch
            disp('Error. Line 552.');
        end
    end
    end
end 

% Find maximum VS across clusters
try
    inc = find(VS_sorted(m,:)<0.99);
    inc2 = find(VS_sorted2(m,:)<0.99);
catch
    disp('Error. Line 562.');
end

try
    VS_max(m) = max(VS_sorted(m,inc));                  % Maximum VS
    VS_max2(m) = max(VS_sorted2(m,inc2));                  % Maximum VS
catch
    disp('Error. Line 569.');
end

try
    maxind = find(VS_sorted(m,:)==VS_max(m));    % Neuron index
    maxind2 = find(VS_sorted2(m,:)==VS_max2(m));    % Neuron index
    VS_max_ind(m) = maxind(1);VS_max_ind2(m) = maxind2(1);
catch
    disp('Error. Line 577.');
end

try
    numspk_max(m) = size(spikes_sorted{m}{VS_max_ind(m)},2);
    numspk_max2(m) = size(spikes_sorted2{m}{VS_max_ind2(m)},2);
    proportion_numspk_max(m) = numspk_max(m)/length(tr{m});
catch
    disp('Error. Line 585.');
end

try
    VSsc_max(m) = max(VSsc_sorted(m,inc));
    maxind = find(VSsc_sorted(m,:)==VSsc_max(m));
    VSsc_max_ind(m) = maxind(1);
catch
    disp('Error. Line 593.');
end

try
    numspk_max_sc(m) = size(spikes_sorted{m}{VSsc_max_ind(m)},2);
    proportion_numspk_max_sc(m) = numspk_max_sc(m)/length(tr{m});
    samepk_scaled(m) = VS_max_ind(m) == VSsc_max_ind(m);
    proportion_numspk_max2(m) = numspk_max2(m)/length(tr{m});
catch
    disp('Error. Line 602.');
end

try
    VSsc_max2(m) = max(VSsc_sorted2(m,inc2));
    maxind2 = find(VSsc_sorted2(m,:)==VSsc_max2(m));
    VSsc_max_ind2(m) = maxind2(1);
catch
    disp('Error. Line 610.');
end

try
    numspk_max_sc2(m) = size(spikes_sorted2{m}{VSsc_max_ind2(m)},2);
    proportion_numspk_max_sc2(m) = numspk_max_sc2(m)/length(tr{m});
    samepk_scaled2(m) = VS_max_ind2(m) == VSsc_max_ind2(m);
catch
    disp('Error. Line 618.');
end


    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Plotting                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear h1;

if plotyn == 1
    disp('Plotting...')
    
    h1 = [];
    colors={'r','g','b','c','o','y','m','k'};

% Plot the raw time series
figure('Name','Time Series');
h1 = [h1 1];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        xc = Xc{nonraw(m)};
        x = Ir{nonraw(m)};
        tvec2 = linspace(0,length(x)/sr,length(x));
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        plot(tvec2,x);hold on;plot(tvec2,xc.*0.01+40,'r')
        xlabel('s');ylabel('pA');title(['Ind=' num2str(m)]);axis off
        end
    end
end

%%%%%%%%%%%
%   PCA   %
%%%%%%%%%%%

% Plot all of the sorted spikes
figure('Name','PCA: Sorted Spikes');
h1 = [h1 2];
for j = 1:length(nonraw)
    try
        if numclust(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean{j});
        for k = 1:numclust(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            plot(tvec{1}(1:size(spikes_sorted{j}{k},1)).*1e3,spikes_sorted{j}{sortB(k)},colors{k});box off;hold all;
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end



% Plot the sorting scatterplots
figure('Name','PCA: Clusters');
h1 = [h1 3];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            scatter(c{m}(T{m}==k,1),c{m}(T{m}==k,2),colors{k});
            if numclust(m) > 1
                title([num2str(numclust(m)) ' clusters'])
            else
                title('1 cluster');
            end
            xlabel(['PC1 (' num2str(explained{m}(1)) '%)']);ylabel(['PC2 (' num2str(explained{m}(2)) '%)']);
            box off;axis off;
            hold all;
        end
        end
    end
end

% Plot the sorted stimuli associated with each spike
figure('Name','PCA: Sorted Stimuli');
h1 = [h1 4];
for j = 1:length(nonraw)
    try
        if numclust(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean{j});
        for k = 1:numclust(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_stim_sorted{j}{k},1)).*1e3,mean(spikes_stim_sorted{j}{sortB(k)},2),std(spikes_stim_sorted{j}{sortB(k)},[],2)./sqrt(size(spikes_stim_sorted{j}{k},1)),colors{k});box off;hold all;
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end

% Plot the mean spike shape for each sorted set
figure('Name','PCA: Mean Sorted Spike Shape');
h1 = [h1 5];
for j = 1:length(nonraw)
    try
        if numclust(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean{j});
        for k = 1:numclust(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_sorted{j}{k},1)).*1e3,mean(spikes_sorted{j}{sortB(k)},2),std(spikes_sorted{j}{sortB(k)},[],2)./sqrt(size(spikes_sorted{j}{k},1)),colors{k});box off;hold all;
            axis([0 4 -max(pkpkampl_mean/1.8) max(pkpkampl_mean/1.8)])
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end

% Plot the mean spike power for each sorted set
figure('Name','PCA: Mean Sorted Spike Power');
h1 = [h1 6];
for j = 1:length(nonraw)
    try
    for k = 1:numclust(j)
        try
            maxx(j,k) = max(mean(spikes_sorted{j}{(k)},2).^2);
        end
    end
    end
end
for j = 1:length(nonraw)
    try
        if numclust(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean{j});
        for k = 1:numclust(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_sorted{j}{k},1)).*1e3,mean(spikes_sorted{j}{sortB(k)},2).^2,std(spikes_sorted{j}{sortB(k)},[],2)./sqrt(size(spikes_sorted{j}{k},1)),colors{k});box off;hold all;
            axis([0 4 0 max(max(maxx))])
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end


% Plot the phase portraits
figure('Name','PCA: Phase Portraits');
h1 = [h1 7];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            p = rose(xc_phase_spikes_sorted{m}{k});            
            t = findall(gca,'type','text');
            delete(t);
            title(['VSmax=' num2str(VS_max(m))])
            box off;axis off;
            hold all;
        end
        end
    end
end

% Plot the ISI distributions
figure('Name','PCA: Interspike Interval Distributions');
h1 = [h1 8];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            [a1, b1] = hist(IEI_sorted{m}{k},freedmandiaconis(IEI_sorted{m}{k}));
            a1 = a1./trapz(b1,a1);
            plot(b1,a1);
            hold all;
        end
        end
    end
end

% Plot the hazard functions
figure('Name','PCA: Hazard Functions');
h1 = [h1 17];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            [a1, b1] = hist(IEI_sorted{m}{k},freedmandiaconis(IEI_sorted{m}{k}));
            a1 = a1./trapz(b1,a1);
            ISI_dist{m}{k}(:,1) = b1;
            ISI_dist{m}{k}(:,2) = a1;
            survivor_func_1{m}{k}(:,1) = b1;
            survivor_func_1{m}{k}(:,2) = 1 - cumtrapz(b1,a1);
            hazard_func_1{m}{k}(:,1) = b1;
            hazard_func_1{m}{k}(:,2) = -a1 ./ survivor_func_1{m}{k}(:,2)';
            plot(hazard_func_1{m}{k}(:,1),hazard_func_1{m}{k}(:,2));
            ylim(0,1);
            hold all;
        end
        end
    end
end

% Plot the survivor functions
figure('Name','PCA: Survivor Functions');
h1 = [h1 18];
for m = 1:length(nonraw)
    try
        if numclust(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            plot(survivor_func_1{m}{k}(:,1),survivor_func_1{m}{k}(:,2));
            ylim([0 1])            
            hold all;
        end
        end
    end
end

%%%%%%%%%%%
% WAVELET %
%%%%%%%%%%%

% Plot all of the sorted spikes (wavelet)
figure('Name','Wavelet: Sorted Spikes');
h1 = [h1 9];
for j = 1:length(nonraw)
    try
        if numclust2(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
        for k = 1:numclust2(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            plot(tvec{1}(1:size(spikes_sorted2{j}{k},1)).*1e3,spikes_sorted2{j}{sortB(k)},colors{k});box off;hold all;
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end


% Plot the sorting scatterplots
figure('Name','Wavelet: Clusters');
h1 = [h1 10];
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust2(m)
            scatter3(cA{m}(coeff1(m),T2{m}==k),cA{m}(coeff2(m),T2{m}==k),cA{m}(coeff3(m),T2{m}==k),colors{k});
            if numclust2(m) > 1
                title([num2str(numclust2(m)) ' clusters'])
            else
                title('1 cluster');
            end
            xlabel('C1');ylabel('C2');
            box off;axis off;
            hold all;
        end
        end
    end
end


% Plot the sorted stimuli associated with each spike
figure('Name','Wavelet: Sorted Stimuli');
h1 = [h1 11];
for j = 1:length(nonraw)
    try
        if numclust2(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
        for k = 1:numclust2(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_stim_sorted2{j}{k},1)).*1e3,mean(spikes_stim_sorted2{j}{sortB(k)},2),std(spikes_stim_sorted2{j}{sortB(k)},[],2)./sqrt(size(spikes_stim_sorted2{j}{k},1)),colors{k});box off;hold all;
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end

% Plot the mean spike shape for each sorted set
figure('Name','Wavelet: Mean Sorted Spike Shape');
h1 = [h1 12];
for j = 1:length(nonraw)
    try
        if numclust2(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
        for k = 1:numclust2(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_sorted2{j}{k},1)).*1e3,mean(spikes_sorted2{j}{sortB(k)},2),std(spikes_sorted2{j}{sortB(k)},[],2)./sqrt(size(spikes_sorted2{j}{k},1)),colors{k});box off;hold all;
            axis([0 4 -max(pkpkampl_mean/1.8) max(pkpkampl_mean/1.8)])
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end

% Plot the mean spike power for each sorted set
figure('Name','Wavelet: Mean Sorted Spike Power');
h1 = [h1 13];
for j = 1:length(nonraw)
    try
    for k = 1:numclust2(j)
        try
            maxx2(j,k) = max(mean(spikes_sorted2{j}{(k)},2).^2);
        end
    end
    end
end
for j = 1:length(nonraw)
    try
        if numclust2(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
        for k = 1:numclust2(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            errorbar(tvec{1}(1:size(spikes_sorted2{j}{k},1)).*1e3,mean(spikes_sorted2{j}{sortB(k)},2).^2,std(spikes_sorted2{j}{sortB(k)},[],2)./sqrt(size(spikes_sorted2{j}{k},1)),colors{k});box off;hold all;
            axis([0 4 0 max(max(maxx2))])
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end


% Plot the phase portraits
figure('Name','Wavelet: Phase Portraits');
h1 = [h1 14];
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust2(m)
            p = rose(xc_phase_spikes_sorted2{m}{k});            
            t = findall(gca,'type','text');
            delete(t);
            title(['VSmax=' num2str(VS_max2(m))])
            box off;axis off;
            hold all;
        end
        end
    end
end

% Plot the ISI distributions
figure('Name','Wavelet: Interspike Interval Distributions');
h1 = [h1 15];
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust2(m)
            [a1, b1] = hist(IEI_sorted2{m}{k},freedmandiaconis(IEI_sorted2{m}{k}));
            a1 = a1./trapz(b1,a1);
            plot(b1,a1);
            hold all;
        end
        end
    end
end


% Plot the hazard functions
figure('Name','Wavelet: Hazard Functions');
h1 = [h1 19];
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust(m)
            [a1, b1] = hist(IEI_sorted2{m}{k},freedmandiaconis(IEI_sorted2{m}{k}));
            a1 = a1./trapz(b1,a1);
            ISI_dist2{m}{k}(:,1) = b1;
            ISI_dist2{m}{k}(:,2) = a1;
            survivor_func_2{m}{k}(:,1) = b1;
            survivor_func_2{m}{k}(:,2) = 1 - cumtrapz(b1,a1);
            hazard_func_2{m}{k}(:,1) = b1;
            hazard_func_2{m}{k}(:,2) = -a1 ./ survivor_func_2{m}{k}(:,2)';
            plot(hazard_func_2{m}{k}(:,1),hazard_func_2{m}{k}(:,2));
            ylim([0 1]);
            hold all;
        end
        end
    end
end

% Plot the survivor functions
figure('Name','Wavelet: Survivor Functions');
h1 = [h1 16];
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust2(m)
            plot(survivor_func_2{m}{k}(:,1),survivor_func_2{m}{k}(:,2));
            ylim([0 1])
            hold all;
        end
        end
    end
end

else
    disp('Skipping plots...');
end

spc1 = 0;

if spc1 == 1
% Supraparamagnetic clustering
disp('Doing supraparamagnetic clustering...');
samplingInterval = 1/Fs;
data = data';
chan = 1;
startData = 0;
par.sr = Fs;sr = Fs;
par.fmin_detect = 0.5*Fs;
par.fmax_detect = 0.1*Fs;
fmin_detect = par.fmin_detect;
fmax_detect = par.fmax_detect;
data = [];
for j = 1:length(Ir)
    data = [data Ir{j}];
end
if datapath(end) == '/'
    save([datapath 'Extracted Data-waveclus.mat']);
    Get_spikes([datapath 'Extracted Data-waveclus.mat']);
    Do_clustering([datapath 'Extracted Data-waveclus_spikes.mat']);
    importdata([datapath 'fig2print_Extracted Data-wave_clus.png']);
    figure;h1 = [h1 17];
    imshow(fig2print_ExtractedData0x2Dwave_clus);
else
    save([datapath '/Extracted Data-waveclus.mat']);
    Get_spikes([datapath '/Extracted Data-waveclus.mat']);
    Do_clustering([datapath '/Extracted Data-waveclus_spikes.mat']);
    importdata([datapath '/fig2print_Extracted Data-wave_clus.png']);
    figure;h1 = [h1 17];
    imshow(fig2print_ExtractedData0x2Dwave_clus);
end
end


% Save the data
if saveyn == 1
    disp('Saving Data...');
    if datapath(end) == '/'
        save([datapath 'Analyzed Data-' num2str(date) '.mat']);
        disp(['Saved as ' datapath 'Analyzed Data-' num2str(date) '.mat'])
        if spc1 == 1
            load([datapath 'times_Extracted Data-wave_clus.mat']);
            importdata([datapath 'spc_log.txt']);loginfo = ans.textdata;
            save([datapath 'Analyzed Data-supraparamagnetic-' num2str(date) '.mat'],'cluster_class','forced','index','inspk','ipermut','par','spikes','Temp','loginfo');
            disp('Saved SPC data.');
        end
    else
        save([datapath '/Analyzed Data-' num2str(date) '.mat']);
        disp(['Saved as ' datapath '/Analyzed Data-' num2str(date) '.mat'])
        if spc1 == 1
            load([datapath '/times_Extracted Data-wave_clus.mat']);
            importdata([datapath '/spc_log.txt']);loginfo = ans.textdata;
            save([datapath '/Analyzed Data-supraparamagnetic-' num2str(date) '.mat'],'cluster_class','forced','index','inspk','ipermut','par','spikes','Temp','loginfo');
            disp('Saved SPC data.');
        end
    end
else
    disp('Data not saved.')
end


% Save the figures
if savefigyn == 1
    disp('Saving Figures...');
    if datapath(end) == '/'
        savefig(h1,[datapath 'Analyzed Data-' num2str(date) '.fig']);
        disp(['Saved as ' datapath 'Analyzed Data-' num2str(date) '.fig'])
    else
        savefig(h1,[datapath '/Analyzed Data-' num2str(date) '.fig']);
        disp(['Saved as ' datapath '/Analyzed Data-' num2str(date) '.fig'])
    end
else
    disp('Figures not saved.')
end



disp('Finished.')

end
