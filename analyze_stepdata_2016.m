% Load data
fn1 = '2016_06_27_0000.abf';

cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\'
cd(fn1(1:end-4));
filenames = dir('Analyzed*.mat');
load(filenames(end).name);

%%
% Extract information
for j = 1:length(Xc)
    disp(num2str(j))
    % Get pulse information
    pulsestart = 1000;
    pulseend = 3000;
    pre = mean(Xc{j}(1:500));
    post = mean(Xc{j}(1500:2000));
    pulse_magnitude(j) = post-pre;
   
    % Get spike rates before, during, and after pulse
    spikerate_prepulse(j) = length(intersect(find(tr{j}>=1), find(tr{j}<=pulsestart-1)))/(tvec{j}(pulsestart-1)-tvec{j}(1));
    spikerate_pulse(j) = length(intersect(find(tr{j}>=pulsestart), find(tr{j}<=pulseend-(pulseend-pulsestart)/2)))/(tvec{j}(pulseend-(pulseend-pulsestart)/2)-tvec{j}(pulsestart));
    spikerate_postpulse(j) = length(intersect(find(tr{j}>=pulseend), find(tr{j}<=length(Xc{j})-(length(Xc{j})-pulseend)/2)))/(tvec{j}(length(Xc{j}-(length(Xc{j})-pulseend)/2))-tvec{j}(pulseend));

    % Onset to first spike
    qr = intersect(find(tr{j}>=pulsestart), find(tr{j}<=pulseend-1));qr=qr(1);
    spikeonset_pulse(j) = tvec{j}(qr) - tvec{j}(pulsestart);
    qr = intersect(find(tr{j}>=pulseend), find(tr{j}<=length(Xc{j})));qr=qr(1);
    spikeonset_postpulse(j) = tvec{j}(qr) - tvec{j}(pulseend);
   
    % Calculate PSTH for spikes
    nn = 1:length(Xc{j});
    qr = findnearest(tvec{j},0.1); qr=qr(1); % 100 ms windows
    rt=buffer(nn,qr,round(qr*0.95));                % 95% overlap
    for k = 1:size(rt,2)
        xc0 = Xc{j}(ceil(mean(rt(:,k))));
        psth_times{j}(k) = tvec{j}(ceil(mean(rt(:,k))));
        psth_spikes{j}(k) = length(intersect(find(tr{j}>=rt(1,k)), find(tr{j}<=rt(end,k))));
    end
   
    nn = 1:pulsestart-1;
    qr = findnearest(tvec{j},0.1); qr=qr(1); % 100 ms windows
    rt=buffer(nn,qr,round(qr*0.95));                % 95% overlap
    for k = 1:size(rt,2)
        xc1 = Xc{j}(ceil(mean(rt(:,k))));
        psth_times_prepulse{j}(k) = tvec{j}(ceil(mean(rt(:,k))));
        psth_spikes_prepulse{j}(k) = length(intersect(find(tr{j}>=rt(1,k)), find(tr{j}<=rt(end,k))));
    end
   
    nn = pulsestart:pulseend-1;
    qr = findnearest(tvec{j},0.1); qr=qr(1); % 100 ms windows
    rt=buffer(nn,qr,round(qr*0.95));                % 95% overlap
    for k = 1:size(rt,2)
        xc2 = Xc{j}(ceil(mean(rt(:,k))));
        psth_times_pulse{j}(k) = tvec{j}(ceil(mean(rt(:,k))));
        psth_spikes_pulse{j}(k) = length(intersect(find(tr{j}>=rt(1,k)), find(tr{j}<=rt(end,k))));
    end
   
    nn = pulseend:length(Xc{j});
    qr = findnearest(tvec{j},0.1); qr=qr(1); % 100 ms windows
    rt=buffer(nn,qr,round(qr*0.95));                % 95% overlap
    for k = 1:size(rt,2)
        xc3 = Xc{j}(ceil(mean(rt(:,k))));
        psth_times_postpulse{j}(k) = tvec{j}(ceil(mean(rt(:,k))));
        psth_spikes_postpulse{j}(k) = length(intersect(find(tr{j}>=rt(1,k)), find(tr{j}<=rt(end,k))));
    end
   
    % Calculate mutual information for entire trace and around different
    % parts of the pulse
    
%     psth_mutualinfo(j) = mutualinformation4(psth_spikes{j},xc0,freedmandiaconis(xc0));
%     psth_mutualinfo_prepulse(j) = mutualinformation4(psth_spikes_prepulse{j},xc1,freedmandiaconis(xc0));
%     psth_mutualinfo_pulse(j) = mutualinformation4(psth_spikes_pulse{j},xc2,freedmandiaconis(xc0));
%     psth_mutualinfo_postpulse(j) = mutualinformation4(psth_spikes_postpulse{j},xc3,freedmandiaconis(xc0));
   
end

%%
pulsemags = unique(pulse_magnitude);

for j = 1:length(pulsemags)
    qr = find(pulse_magnitude==pulsemags(j));
    qrn = find(pulse_magnitude==-pulsemags(j));
    
    meanspikerate(j) = mean(horzcat(spikerate_pulse(qr),spikerate_postpulse(qrn)));
    semspikerate(j) = mean(horzcat(spikerate_pulse(qr),spikerate_postpulse(qrn)))./sqrt(length(qr)+length(qrn));
    meanspikeonset(j) = mean(horzcat(spikeonset_pulse(qr),spikeonset_postpulse(qrn)));
    semspikeonset(j) = mean(horzcat(spikeonset_pulse(qr),spikeonset_postpulse(qrn)))./sqrt(length(qr)+length(qrn));
    
    for k = 1:length(qr)
        rasterplot_x{j}{k} = tvec{1}(tr{qr(k)});
        rasterplot_y{j}{k} = k.*ones(1,length(tr{qr(k)}));
    end
end
    
   
