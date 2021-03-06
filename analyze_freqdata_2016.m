fL = length(Ir)/length(freqstim)
for j = 1:length(freqstim)
    for k = 1:fL
        VS1(j,k) = VS(j*k);
        VS1_max(j,k) = max(VS_sorted(j*k,:));
        VS1_max2(j,k) =max(VS_sorted2(j*k,:));
        qr = find(VS_sorted(j*k,:)==max(VS_sorted(j*k,:)));
        qr2 = find(VS_sorted2(j*k,:)==max(VS_sorted2(j*k,:)));
        mean_phase1(j,k) = meanphase(j*k);
        try
        mean_phase1_sorted(j,k) = meanphase_sorted{j*k}(qr);
        catch
            mean_phase1_sorted(j,k) = NaN;
        end
        try
        mean_phase1_sorted2(j,k) = meanphase_sorted2{j*k}(qr2);
        catch
            mean_phase1_sorted2(j,k) = NaN;
        end
    end
    
    VSmean(j) = nanmean(VS1(j,:));
    VSsem(j) = nanstd(VS1(j,:))/sqrt(fL);
    VSmean_sorted(j) = nanmean(VS1_max(j,:));
    VSsem_sorted(j) = nanstd(VS1_max(j,:))/sqrt(fL);
    VSmean_sorted2(j) = nanmean(VS1_max2(j,:));
    VSsem_sorted2(j) = nanstd(VS1_max2(j,:))/sqrt(fL);
    mean_phase_mean(j) = nanmean(mean_phase1(j,:));
    mean_phase_sem(j) = nanstd(mean_phase1(j,:))/sqrt(fL);
    mean_phase_sorted_mean(j) = nanmean(mean_phase1_sorted(j,:));
    mean_phase_sorted_sem(j) = nanstd(mean_phase1_sorted(j,:))/sqrt(fL);
    mean_phase_sorted2_mean(j) = nanmean(mean_phase1_sorted2(j,:));
    mean_phase_sorted2_sem(j) = nanstd(mean_phase1_sorted2(j,:))/sqrt(fL);
    
end

figure;
subplot(2,3,1);errorbar(freqstim,VSmean,VSsem);ylabel('VS');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean)*0.5 max(VSmean)*1.5]);
subplot(2,3,2);errorbar(freqstim,VSmean_sorted,VSsem_sorted);ylabel('VS');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean_sorted)*0.5 max(VSmean_sorted)*1.5]);
subplot(2,3,3);errorbar(freqstim,VSmean_sorted2,VSsem_sorted2);ylabel('VS');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean_sorted2)*0.5 max(VSmean_sorted2)*1.5]);
subplot(2,3,4);errorbar(freqstim,mean_phase_mean,mean_phase_sem);ylabel('phase');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_mean)*0.5 max(mean_phase_mean)*1.5]);
subplot(2,3,5);errorbar(freqstim,mean_phase_sorted_mean,mean_phase_sorted_sem);ylabel('phase');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_sorted_mean)*0.5 max(mean_phase_sorted_mean)*1.5]);
subplot(2,3,6);errorbar(freqstim,mean_phase_sorted2_mean,mean_phase_sorted2_sem);ylabel('phase');xlabel('Stimulus frequency (Hz)')
% axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_sorted2_mean)*0.5 max(mean_phase_sorted_mean)*1.5]);

%% 
setfiguredefaults(5);

inds = 1:length(freqstim);
inds = 1:length(freqstim);
inds = [1:7 9 11:length(freqstim)];

% close all;


figure
h=errorbar(freqstim(inds),VSmean_sorted2(inds),VSsem_sorted2(inds));
xlabel('Stimulus frequency');
ylabel('Vector strength')
ha = get(h,'Parent');
set(ha,'XScale','log')
set(ha,'XLim',[min(freqstim(inds))-0.2*min(freqstim(inds)) max(freqstim(inds))+0.2*max(freqstim(inds))])

