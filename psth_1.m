
% Calculate peri-stimulus time histograms for complete time series,
% including those from spike sorting by Wavelet and PCA

perc_overlap = 0.9;         % percent overlap
winlength = 0.05;           % window length (sec)
winlength = winlength*sr;   % window length (points)

for m = 1:length(nonraw)
    try
    clear xrt
    disp(num2str(m));
    
    xc = Xc{nonraw(m)}-smooth(Xc{nonraw(m)},2e4);
    xd = Xd{nonraw(m)}-smooth(Xd{nonraw(m)},2e4);
    
    overlap(m) = winlength - winlength*perc_overlap;
    numwin_psth(m) = (length(xc) - winlength) / overlap(m) + 1;
    [xcsplit,inds,indf] = enframe(xc,winlength,overlap(m));
    inc = overlap(m);
    nx = length(xc);
    nf = fix((nx-winlength+inc)/inc);
    xrt = indf(:,ones(1,winlength))+inds(ones(nf,1),:);
    sizex = size(xrt);
    for j = 1:sizex(1)
%         psth{m}(j) = length(intersect(find(pk{m}>1+(j-1)*length(xc)/numwin_psth(m)),find(pk{m}<(j)*length(xc)/numwin_psth(m))));
        psth{m}(j) = length(intersect(find(pk{m} > xrt(j,1)),find(pk{m} < xrt(j,sizex(2)))));
    end
    for k = 1:numclust(m)
        for j = 1:numwin_psth(m)
%             psth_sorted{m,k}(j) = length(intersect(find(tr_sorted{m}{k}>1+(j-1)*length(xc)/numwin_psth(m)),find(tr_sorted{m}{k}<(j)*length(xc)/numwin_psth(m))));
            psth_sorted{m,k}(j) = length(intersect(find(tr_sorted{m}{k} > xrt(j,1)),find(tr_sorted{m}{k} < xrt(j,sizex(2)))));
        end
        
        [pxx_psth_sorted{m,k},fxx_psth_sorted{m,k}]=pwelch(psth_sorted{m,k},[],[],[],sr);
        [pxx_stim_sorted{m,k},fxx_stim_sorted{m,k}]=pwelch(xc,[],[],[],sr);
        qr = find(pxx_stim_sorted{m,k} == max(pxx_stim_sorted{m,k}));qr=qr(1);
        stim_freq_pxx_sorted(m,k) = fxx_stim_sorted{m,k}(qr);
        stim_ampl_pxx_sorted(m,k) = pxx_stim_sorted{m,k}(qr);
    
        qr = findnearest(fxx_psth_sorted{m,k}, stim_freq_pxx_sorted(m,k));qr=qr(1);
        psth_ampl_pxx_sorted(m,k) = pxx_psth_sorted{m,k}(qr);
        psth_freq_pxx_sorted(m,k) = fxx_psth_sorted{m,k}(qr);
    end
    psth_ampl_pxx_max_sorted(m) = max(psth_ampl_pxx_sorted(m,:));
    for k = 1:numclust2(m)
        for j = 1:numwin_psth(m)
%             psth_sorted2{m,k}(j) = length(intersect(find(tr_sorted2{m}{k}>1+(j-1)*length(xc)/numwin_psth(m)),find(tr_sorted2{m}{k}<(j)*length(xc)/numwin_psth(m))));
            psth_sorted2{m,k}(j) = length(intersect(find(tr_sorted2{m}{k} > xrt(j,1)),find(tr_sorted2{m}{k} < xrt(j,sizex(2)))));
        end
        [pxx_psth_sorted2{m,k},fxx_psth_sorted2{m,k}]=pwelch(psth_sorted2{m,k},[],[],[],sr/(length(xc)/numwin_psth(m)));
        [pxx_stim_sorted2{m,k},fxx_stim_sorted2{m,k}]=pwelch(xc,[],[],[],sr);
        qr = find(pxx_stim_sorted2{m,k} == max(pxx_stim_sorted2{m,k}));qr=qr(1);
        stim_freq_pxx_sorted2(m,k) = fxx_stim_sorted2{m,k}(qr);
        stim_ampl_pxx_sorted2(m,k) = pxx_stim_sorted2{m,k}(qr);
    
        qr = findnearest(fxx_psth_sorted2{m,k}, stim_freq_pxx_sorted2(m,k));qr=qr(1);
        psth_ampl_pxx_sorted2(m,k) = pxx_psth_sorted2{m,k}(qr);
        psth_freq_pxx_sorted2(m,k) = fxx_psth_sorted2{m,k}(qr);
    end
    psth_ampl_pxx_max_sorted2(m) = max(psth_ampl_pxx_sorted2(m,:));
    
    [pxx_psth{m},fxx_psth{m}]=pwelch(psth{m},[],[],[],sr/(length(xc)/numwin_psth(m)));
    [pxx_stim{m},fxx_stim{m}]=pwelch(xc,[],[],[],sr);
    [pxx_resp{m},fxx_resp{m}]=pwelch(xd,[],[],[],sr);
    
    qr = find(pxx_stim{m} == max(pxx_stim{m}));qr=qr(1);
    stim_freq_pxx(m) = fxx_stim{m}(qr);
    stim_ampl_pxx(m) = pxx_stim{m}(qr);
    
    qr = find(pxx_resp{m} == max(pxx_resp{m}));qr=qr(1);
    resp_freq_pxx(m) = fxx_resp{m}(qr);
    resp_ampl_pxx(m) = pxx_resp{m}(qr);
    
    qr = findnearest(fxx_psth{m}, stim_freq_pxx(m));qr=qr(1);
    psth_ampl_pxx(m) = pxx_psth{m}(qr);
    psth_freq_pxx(m) = fxx_psth{m}(qr);
    
    end
    
end


% Plot Fourier transforms of peri-stimulus time histograms
figure('Name','PSTH FFT');
for j = 1:length(nonraw)
    try
            subplot_tight(ceil(sqrt(length(nonraw))),ceil(sqrt(length(nonraw))),j,[0.04 0.01]);
            loglog(fxx_stim{j},pxx_stim{j},'k');hold on;axis auto
            loglog(fxx_psth{j},pxx_psth{j},'r');xlim([0.1 500]);
    end
end
    
    
