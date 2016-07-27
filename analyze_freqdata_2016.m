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
    
    VSmean(j) = mean(VS1(j,:));
    VSsem(j) = std(VS1(j,:))/sqrt(fL);
    VSmean_sorted(j) = mean(VS1_max(j,:));
    VSsem_sorted(j) = std(VS1_max(j,:))/sqrt(fL);
    VSmean_sorted2(j) = mean(VS1_max2(j,:));
    VSsem_sorted2(j) = std(VS1_max2(j,:))/sqrt(fL);
    mean_phase_mean(j) = mean(mean_phase1(j,:));
    mean_phase_sem(j) = std(mean_phase1(j,:))/sqrt(fL);
    mean_phase_sorted_mean(j) = mean(mean_phase1_sorted(j,:));
    mean_phase_sorted_sem(j) = std(mean_phase1_sorted(j,:))/sqrt(fL);
    mean_phase_sorted2_mean(j) = mean(mean_phase1_sorted2(j,:));
    mean_phase_sorted2_sem(j) = std(mean_phase1_sorted2(j,:))/sqrt(fL);
    
end

figure;
subplot(2,3,1);errorbar(freqstim,VSmean,VSsem);ylabel('VS');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean)*0.5 max(VSmean)*1.5]);
subplot(2,3,2);errorbar(freqstim,VSmean_sorted,VSsem_sorted);ylabel('VS');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean_sorted)*0.5 max(VSmean_sorted)*1.5]);
subplot(2,3,3);errorbar(freqstim,VSmean_sorted2,VSsem_sorted2);ylabel('VS');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(VSmean_sorted2)*0.5 max(VSmean_sorted2)*1.5]);
subplot(2,3,4);errorbar(freqstim,mean_phase_mean,mean_phase_sem);ylabel('phase');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_mean)*0.5 max(mean_phase_mean)*1.5]);
subplot(2,3,5);errorbar(freqstim,mean_phase_sorted_mean,mean_phase_sorted_sem);ylabel('phase');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_sorted_mean)*0.5 max(mean_phase_sorted_mean)*1.5]);
subplot(2,3,6);errorbar(freqstim,mean_phase_sorted2_mean,mean_phase_sorted2_sem);ylabel('phase');
axis([min(freqstim)*0.5 max(freqstim)*1.5 mean(mean_phase_sorted2_mean)*0.5 max(mean_phase_sorted_mean)*1.5]);

