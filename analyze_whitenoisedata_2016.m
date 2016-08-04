clear spikes_stim_sorted_mean
clear spikes_stim_sorted_sem
clear spikes_stim_sorted2_mean
clear spikes_stim_sorted2_sem

for j = 1:length(nonraw)
[~,sortB] = sort(spikes_sorted_pkpkampl_mean{j});
for k = 1:numclust(j)
spikes_stim_sorted_mean{j}(:,k) = mean(spikes_stim_sorted{j}{sortB(k)},2);
spikes_stim_sorted_sem{j}(:,k) = std(spikes_stim_sorted{j}{sortB(k)},[],2)./sqrt(size(spikes_stim_sorted{j}{k},1));
end
end
for j = 1:length(nonraw)
[~,sortB] = sort(spikes_sorted_pkpkampl_mean2{j});
for k = 1:numclust2(j)
spikes_stim_sorted2_mean{j}(:,k) = mean(spikes_stim_sorted2{j}{sortB(k)},2);
spikes_stim_sorted2_sem{j}(:,k) = std(spikes_stim_sorted2{j}{sortB(k)},[],2)./sqrt(size(spikes_stim_sorted2{j}{k},1));
end
end

for j = 1:length(nonraw)
    try
    qr(j) = find(VS_sorted(j,:)==VS_max(j));
    catch
        qr(j)=1;
    end
    try
    qr2(j) = find(VS_sorted2(j,:)==VS_max2(j));
    catch
        qr2(j)=1;
    end
%     vars1 = var(spikes_stim_sorted_mean{j}(1:200,:));qr(j) = find(vars1==max(vars1));
%     vars2 = var(spikes_stim_sorted2_mean{j}(1:200,:));qr2(j) = find(vars2==max(vars2));
    try
    stimulus_mean(:,j) = spikes_stim_sorted_mean{j}(:,qr(j));
    end
    try
    stimulus_mean2(:,j) = spikes_stim_sorted2_mean{j}(:,qr2(j));
    end
    
end
    
%%
LL = size(stimulus_mean,2);
figure;
for j = 1:LL
    subplot_tight(ceil(sqrt(LL)),ceil(sqrt(LL)),j);
    plot((dt-length(stimulus_mean2)/2*dt:dt:length(stimulus_mean2)/2*dt).*1e3,mean(stimulus_mean2(:,j),2),'r');
    plot((dt-length(stimulus_mean)/2*dt:dt:length(stimulus_mean)/2*dt).*1e3,mean(stimulus_mean(:,j),2),'k');
    try
    axis([0 length(stimulus_mean)/2*dt.*1e3 -max(max(stimulus_mean(:,j))) max(max(stimulus_mean(:,j)))]);
    end
end
%%
inds = 1:2:LL;
figure;
plot((dt-length(stimulus_mean2)/2*dt:dt:length(stimulus_mean2)/2*dt).*1e3,mean(stimulus_mean2(:,inds),2))

cov1= cov([mean(stimulus_mean2(1:250,6),2)';mean(stimulus_mean2(1:250,7),2)';mean(stimulus_mean2(1:250,8),2)']);
cov2= cov([mean(stimulus_mean2(251:500,6),2)';mean(stimulus_mean2(251:500,7),2)';mean(stimulus_mean2(251:500,8),2)']);