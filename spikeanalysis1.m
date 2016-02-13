%% Load data
clear all;close all;clc
load('/Users/joshsalvi/Documents/Lab/Lab/Ephys Data/2016-02-05.01/Ear 1/Cell 1/Extracted Data.mat')

%% Spike analysis


manualyn = 0;
setfiguredefaults;

% initialization
for m = 1:length(nonraw)
    disp(['Iteration ' num2str(m) ' out of ' num2str(length(nonraw))]);
    
        % Define command signal and neural response
        xc = Xc{nonraw(m)};
        x = Ir{nonraw(m)};
        
        % Sample rate (Hz)
        sr = 10e3;
        
        % Spike length for extraction (s)
        sl = 0.004;
        sl = sl*sr;
        
        % Extract spikes, peaks in command signal, set maximum frequency
        % limit for optimization purposes
        thr = 10*median(abs(x)./0.6745);
        [pk{m}, tr{m}] = PTDetect(x, thr);
        [pkstim{m}, trstim{m}] = PTDetect(xc, std(xc-mean(xc)));
        numpkstim = length(pkstim{m});
        sl_stim = length(xc)/numpkstim;
        stim_freq(m) = numpkstim/length(xc)*10e3;
        max_freq = 1e3;
        

if length(tr{m}) > 0 && length(tr{m}) < 5e3
    for j = 1:length(tr{m})  
        
            % Extract spikes
            try
                spikes{m}(:,j) = x(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2));
                spikes_stim{m}(:,j) = xc(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2)); 
            catch
                spikes{m}(:,j) = zeros(1,length(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2)));
                spikes_stim{m}(:,j) = zeros(1,length(tr{m}(j) - round(sl/2):tr{m}(j) + round(sl/2)));                                 
            end      
            
            % Extract spike amplitudes, power, and areas
            pkpkampl{m}(j) = max(spikes{m}(:,j)) - min(spikes{m}(:,j));
            spike_power{m}(j) = trapz(tvec{m}(1:length(spikes{m}(:,j))),spikes{m}(:,j).^2);
            spike_area{m}(j) = trapz(tvec{m}(1:length(spikes{m}(:,j))),spikes{m}(:,j));
                  
    end
    
    % Mean amplitude (all spikes)
    pkpkampl_mean(m) = mean(pkpkampl{m});
    pkpkampl_sem(m) = std(pkpkampl{m})/sqrt(length(tr{m}));
    
    % Mean power (all spikes)
    spike_power_mean(m) = mean(spike_power{m});
    spike_power_sem(m) = std(spike_power{m})/sqrt(length(tr{m}));
    
    % Mean area (all spikes)
    spike_area_mean(m) = mean(spike_area{m});
    spike_area_sem(m) = std(spike_area{m})/sqrt(length(tr{m}));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting by Principal components analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = 1;
    clear nc dim1
    
    try
        nc = numclust(m);
    end
    try 
        dim1 = dims{m};
    end
    m2=1;
    if exist('dim1') == 0
        dims{m} = 2;
    end
     while r == 1 && exist('nc') == 0 
         %}
        
        [c{m}, score{m}, latent{m}, tsquared{m}, explained{m}, mu{m}]=pca(spikes{m});
        
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r');
        subplot(2,1,2);scatter(c{m}(:,1),c{m}(:,2));
        xlabel(['PC1 (' num2str(explained{m}(1)) '%)']);ylabel(['PC2 (' num2str(explained{m}(2)) '%)']);
        
        
        clear clust1
        
        if manualyn == 1
            if m2 > 1
                numclust(m) = input('Number of clusters:   ');
                dims{m} = input('Dimensions (1,2,3,...):   ');
            else
                for j = 2:5
%                     clust1(:,j-1) = kmeans(c{m}(:,1:dims{m}),j,'Replicates',10);
                end
                evalclust{m} = evalclusters(c{m}(:,1:dims{m}),'kmeans','gap','KList',1:4);
                numclust(m) = evalclust{m}.OptimalK;
            end
        else
            for j = 2:5
%                 clust1(:,j-1) = kmeans(c{m}(:,1:dims{m}),j,'Replicates',10);
            end
            evalclust{m} = evalclusters(c{m}(:,1:dims{m}),'kmeans','gap','KList',1:4);
            numclust(m) = evalclust{m}.OptimalK;
        end
        
        T{m} = kmeans(c{m}(:,1:dims{m}),numclust(m),'Replicates',10);
        explainedvariance(m) = sum(explained{m}(1:dims{m}));
        
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r')
        
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
        if manualyn == 1
            r = input('Repeat? (1=yes):   ');
        else
            r = 0;
        end
        m2=m2+1;
    
    end
    
    
    try
        [c{m}, score{m}, latent{m}, tsquared{m}, explained{m}, mu{m}]=pca(spikes{m});
        T{m} = kmeans(c{m}(:,1:dims{m}),numclust(m),'Replicates',10);
    end
%}
    
    % Sort the spikes and calculate statistics for each subset
    if numclust(m) >=2
    for k = 1:numclust(m)
        spikes_sorted{m}{k} = spikes{m}(:,T{m}==k);
        spikes_stim_sorted{m}{k} = spikes_stim{m}(:,T{m}==k);
        tr_sorted{m}{k} = tr{m}(T{m}==k);
%         pk_sorted{m}{k} = pk{m}(T{m}==k);
    end
    else
        spikes_sorted{m}{1} = spikes{m};
        spikes_stim_sorted{m}{k} = spikes_stim{m};
        tr_sorted{m}{1} = tr{m};
%         pk_sorted{m}{1} = pk{m};
    end
    
    for k = 1:numclust(m)
        spikes_sorted_pow_mean{m}(k) = mean(spike_power{m}(T{m}==k));
        spikes_sorted_pow_sem{m}(k) = std(spike_power{m}(T{m}==k))/sqrt(sum(T{m}==k));
        spikes_sorted_area_mean{m}(k) = mean(spike_area{m}(T{m}==k));
        spikes_sorted_area_sem{m}(k) = std(spike_area{m}(T{m}==k))/sqrt(sum(T{m}==k));
        spikes_sorted_pkpkampl_mean{m}(k) = mean(pkpkampl{m}(T{m}==k));
        spikes_sorted_pkpkampl_sem{m}(k) = std(pkpkampl{m}(T{m}==k))/sqrt(sum(T{m}==k));
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort spikes using the WAVELET transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:size(spikes{m},2)
    % Discrete time wavelet transform
    [cA{m}(:,k),cD{m}(:,k)] = dwt(spikes{m}(:,k),'haar');
end
for l = 1:size(cA{m},1)
    % Find the coefficients that show greatest spread using Lilliefors test
    % A large statistic corresponds to a divergence from a normal
    % distribution
    warning off
    [~,~,kstat{m}(l)] = lillietest(cA{m}(l,:));
end

% Sort by K-statistic
sortedK = sort(kstat{m});
coeff1 = find(kstat{m} == sortedK(end));
coeff2 = find(kstat{m} == sortedK(end-1));

% Cluster data
clear nc
r = 1;
    try
        nc = numclust2(m);
    end
    try 
        dim1 = dims{m};
    end
m2=1;
    if exist('dim1') == 0
        dims{m} = 2;
    end
     while r == 1 && exist('nc') == 0 
         %}
        
       
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r');
        subplot(2,1,2);scatter(cA{m}(coeff1,:),cA{m}(coeff2,:));
        xlabel('C1');ylabel('C2');
        
                
        clear clust1 clust2
        
        if manualyn == 1
            if m2 > 1
                numclust2(m) = input('Number of clusters:   ');
            else
                for j = 2:5
%                     clust2(:,j-1) = kmeans([cA{m}(coeff1,:);cA{m}(coeff2,:)]',j,'Replicates',10);
                end
                evalclust2{m} = evalclusters([cA{m}(coeff1,:);cA{m}(coeff2,:)]','kmeans','gap','KList',1:4);
                numclust2(m) = evalclust2{m}.OptimalK;
            end
        else
            for j = 2:5
%                 clust2(:,j-1) = kmeans([cA{m}(coeff1,:);cA{m}(coeff2,:)]',j,'Replicates',10);
            end
            evalclust2{m} = evalclusters([cA{m}(coeff1,:);cA{m}(coeff2,:)]','kmeans','gap','KList',1:4);
            numclust2(m) = evalclust2{m}.OptimalK;
        end
        
        T2{m} = kmeans([cA{m}(coeff1,:);cA{m}(coeff2,:)]',numclust2(m),'Replicates',10);
        
        close all;figure;
        subplot(2,1,1);plot(x);hold on;plot(xc.*0.01+40,'r')
        
        if numclust2(m) >= 2
        colors = {'r','g','b','y','k','m','c'};    
        for k = 1:numclust2(m)
            subplot(2,1,2);scatter(cA{m}(coeff1,T2{m}==k),cA{m}(coeff2,T2{m}==k),colors{k});
            title([num2str(numclust2(m)) ' clusters'])
            xlabel('C1');ylabel('C2');
            hold all;
        end
        
        else
            
            k=1;subplot(2,1,2);scatter(cA{m}(coeff1,:),cA{m}(coeff2,:));
            
        end
        if manualyn == 1
            r = input('Repeat? (1=yes):   ');
        else 
            r=0;
        end
        m2=m2+1;
    
    end
    
    
    try
        T2{m} = kmeans([cA{m}(coeff1,:);cA{m}(coeff2,:)]',numclust2(m),'Replicates',10);
    end
    
% Sort the spikes by clustering from wavelet coefficients
    if numclust2(m) >=2
    for k = 1:numclust2(m)
        spikes_sorted2{m}{k} = spikes{m}(:,T2{m}==k);
        spikes_stim_sorted2{m}{k} = spikes_stim{m}(:,T2{m}==k);
        tr_sorted2{m}{k} = tr{m}(T2{m}==k);
%         pk_sorted{m}{k} = pk{m}(T{m}==k);
    end
    else
        spikes_sorted2{m}{1} = spikes{m};
        spikes_stim_sorted2{m}{k} = spikes_stim{m};
        tr_sorted2{m}{1} = tr{m};
%         pk_sorted{m}{1} = pk{m};
    end
    
    for k = 1:numclust2(m)
        spikes_sorted_pow_mean2{m}(k) = mean(spike_power{m}(T2{m}==k));
        spikes_sorted_pow_sem2{m}(k) = std(spike_power{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
        spikes_sorted_area_mean2{m}(k) = mean(spike_area{m}(T2{m}==k));
        spikes_sorted_area_sem2{m}(k) = std(spike_area{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
        spikes_sorted_pkpkampl_mean2{m}(k) = mean(pkpkampl{m}(T2{m}==k));
        spikes_sorted_pkpkampl_sem2{m}(k) = std(pkpkampl{m}(T2{m}==k))/sqrt(sum(T2{m}==k));
    end
    
% Extract phases

    % Hilbert transform of command signal
    xc_hilb = hilbert(xc);
    xc_phase = atan2(imag(xc_hilb),real(xc_hilb));
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
    
    clear pks_present
    
    % Bernoulli trials: Determine whether a spike occurs in each cycle of
    % the stimulus.
    if stim_freq(m) < max_freq
    for j = 1:numpkstim
        pks_present(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}))>0;
    end
    
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
    end

    
    % Repeat for sorted spikes (PCA)
    for k = 1:numclust(m)
        try
        xc_phase_spikes_sorted{m}{k} = xc_phase_spikes{m}(T{m}==k);
        meanphase_sorted{m}(k) = mean(xc_phase_spikes_sorted{m}{k});
        semphase_sorted{m}(k) = std(xc_phase_spikes_sorted{m}{k})/sqrt(sum(T{m}==k));
        VS_sorted(m,k) = 1/length(tr{m}(T{m}==k))*abs(sum(exp(sqrt(-1).*xc_phase_spikes_sorted{m}{k})));
        if stim_freq(m) < max_freq
        clear pks_present
        for j = 1:numpkstim
            pks_present(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}(T{m}==k)))>0;
        end
        scaling_fac_VS_sorted(m,k) = sum(pks_present)/numpkstim;
        if scaling_fac_VS_sorted(m,k) > 1
            warning(['Scaling factor for vector strength exceeds 1 at index ' num2str(m) '.']);
        elseif scaling_fac_VS_sorted(m,k) < 0
            warning(['Scaling factor for vector strength below 0 at index ' num2str(m) '.']);
        else
            VSsc_sorted(m,k) = VS_sorted(m,k) * scaling_fac_VS_sorted(m,k);
        end
        end
        end
    end
    
    % Repeat for sorted spikes (WAVELET)
    for k = 1:numclust2(m)
        try
        xc_phase_spikes_sorted2{m}{k} = xc_phase_spikes{m}(T2{m}==k);
        meanphase_sorted2{m}(k) = mean(xc_phase_spikes_sorted2{m}{k});
        semphase_sorted2{m}(k) = std(xc_phase_spikes_sorted2{m}{k})/sqrt(sum(T2{m}==k));
        VS_sorted2(m,k) = 1/length(tr{m}(T2{m}==k))*abs(sum(exp(sqrt(-1).*xc_phase_spikes_sorted2{m}{k})));
        if stim_freq(m) < max_freq
        clear pks_present
        for j = 1:numpkstim
            pks_present2(j) = length(intersect(1 + (j-1)*sl_stim:j*sl_stim,tr{m}(T2{m}==k)))>0;
        end
        scaling_fac_VS_sorted2(m,k) = sum(pks_present2(j))/numpkstim;
        if scaling_fac_VS_sorted2(m,k) > 1
            warning(['Scaling factor for vector strength exceeds 1 at index ' num2str(m) '.']);
        elseif scaling_fac_VS_sorted(m,k) < 0
            warning(['Scaling factor for vector strength below 0 at index ' num2str(m) '.']);
        else
            VSsc_sorted2(m,k) = VS_sorted2(m,k) * scaling_fac_VS_sorted2(m,k);
        end
        end
        end
    end
end 

% Find maximum VS across clusters
try
    inc = find(VS_sorted(m,:)<0.99);
    inc2 = find(VS_sorted2(m,:)<0.99);
    VS_max(m) = max(VS_sorted(m,inc));                  % Maximum VS
    VS_max2(m) = max(VS_sorted2(m,inc2));                  % Maximum VS
    maxind = find(VS_sorted(m,:)==VS_max(m));    % Neuron index
    maxind2 = find(VS_sorted2(m,:)==VS_max2(m));    % Neuron index
    VS_max_ind(m) = maxind(1);VS_max_ind2(m) = maxind2(1);
    numspk_max(m) = size(spikes_sorted{m}{VS_max_ind(m)},2);
    numspk_max2(m) = size(spikes_sorted2{m}{VS_max_ind2(m)},2);
    proportion_numspk_max(m) = numspk_max(m)/length(tr{m});
    VSsc_max(m) = max(VSsc_sorted(m,inc));
    maxind = find(VSsc_sorted(m,:)==VSsc_max(m));
    VSsc_max_ind(m) = maxind(1);
    numspk_max_sc(m) = size(spikes_sorted{m}{VSsc_max_ind(m)},2);
    proportion_numspk_max_sc(m) = numspk_max_sc(m)/length(tr{m});
    samepk_scaled(m) = VS_max_ind(m) == VSsc_max_ind(m);
    proportion_numspk_max2(m) = numspk_max2(m)/length(tr{m});
    VSsc_max2(m) = max(VSsc_sorted2(m,inc2));
    maxind2 = find(VSsc_sorted2(m,:)==VSsc_max2(m));
    VSsc_max_ind2(m) = maxind2(1);
    numspk_max_sc2(m) = size(spikes_sorted2{m}{VSsc_max_ind2(m)},2);
    proportion_numspk_max_sc2(m) = numspk_max_sc2(m)/length(tr{m});
    samepk_scaled2(m) = VS_max_ind2(m) == VSsc_max_ind2(m);
end

end

%% PLOTTING
close all;
colors={'r','g','b','c','o','y','m','k'};

% Plot the raw time series
figure('Name','Time Series');
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
set(0,'defaultAxesColorOrder',jet(5))
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


%%%%%%%%%%%
% WAVELET %
%%%%%%%%%%%

% Plot all of the sorted spikes (wavelet)
figure('Name','Wavelet: Sorted Spikes');
for j = 1:length(nonraw)
    try
        if numclust2(j) > 0
        [~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
        for k = 1:numclust(j)
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            plot(tvec{1}(1:size(spikes_sorted2{j}{k},1)).*1e3,spikes_sorted2{j}{sortB(k)},colors{k});box off;hold all;
            xlabel('ms');ylabel('pA');title(['Ind=' num2str(j)]);axis off
        end
        end
    end
end


% Plot the sorting scatterplots
figure('Name','Wavelet: Clusters');
for m = 1:length(nonraw)
    try
        if numclust2(m) > 0
        subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),m,[0.04 0.01]);
        for k = 1:numclust2(m)
            scatter(cA{m}(coeff1,T2{m}==k),cA{m}(coeff2,T2{m}==k),colors{k});
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
set(0,'defaultAxesColorOrder',jet(5))
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
