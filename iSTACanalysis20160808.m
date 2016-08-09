bw_time = 0.001; %sec
ts_time = 0.05;  %sec
bw = Fs*bw_time;    %points
ts = Fs*ts_time;    %points

clear Stim spn sp

disp('extracting spikes...')
sp = zeros(1,length(Xc{1}));
for j = 1:length(tr)
    sp(tr{j}) = sp(tr{j})+1;
end
Stim = Xc{1};
disp('STC...')
[sta,stc,rawmu,rawcov] = simpleSTC(Stim, sp', ts);
[u,s,v] = svd(stc);

ndims = 5;
disp('filters...');
eigvalthresh = 0.05;
[vecs, vals, DD] = compiSTAC(sta, stc, rawmu, rawcov, ndims,eigvalthresh);
KLcontributed = [vals(1); diff(vals)];
ndims = length(vals);


sw = 10;
for j = 1:ndims
    vecs_smoothed(:,j) = smooth(vecs(:,j),sw,'moving');
end

close all;

figure;
subplot(221);
plot(tvec{1}(1:length(sta)).*1e3, sta./norm(sta));
title('STA');xlabel('time (ms)');
axis([0 tvec{1}(length(sta)).*1e3 -max(abs(sta./norm(sta))) max(abs(sta./norm(sta)))]);

subplot(222);  plot(1:ndims, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(223);
plot(tvec{1}(1:length(vecs)).*1e3, vecs(:,1:2));
title('iSTAC filters');xlabel('time (ms)');
legend('1st', '2nd', 'location', 'northwest');
axis([0 tvec{1}(length(sta)).*1e3 -max(max(abs(vecs(:,1:2)))) max(max(abs(vecs(:,1:2))))]);

subplot(224);
plot(tvec{1}(1:length(vecs)).*1e3, vecs(:,3:4));
title('iSTAC filters');xlabel('time (ms)');
legend('3rd', '4th', 'location', 'northwest');
axis([0 tvec{1}(length(sta)).*1e3 -max(max(abs(vecs(:,3:4)))) max(max(abs(vecs(:,3:4))))]);


figure;
subplot(221);
plot(tvec{1}(1:length(sta)).*1e3, sta./norm(sta));
title('STA');xlabel('time (ms)');
axis([0 tvec{1}(length(sta)).*1e3 -max(abs(sta./norm(sta))) max(abs(sta./norm(sta)))]);

subplot(222);  plot(1:ndims, KLcontributed, 'o');
title('KL contribution');
xlabel('subspace dimensionality');

subplot(223);
plot(tvec{1}(1:length(vecs_smoothed)).*1e3, vecs_smoothed(:,1:2));
title('iSTAC filters');xlabel('time (ms)');
legend('1st', '2nd', 'location', 'northwest');
axis([0 tvec{1}(length(sta)).*1e3 -max(max(abs(vecs_smoothed(:,1:2)))) max(max(abs(vecs_smoothed(:,1:2))))]);

subplot(224);
plot(tvec{1}(1:length(vecs_smoothed)).*1e3, vecs_smoothed(:,3:4));
title('iSTAC filters');xlabel('time (ms)');
legend('3rd', '4th', 'location', 'northwest');
axis([0 tvec{1}(length(sta)).*1e3 -max(max(abs(vecs_smoothed(:,3:4)))) max(max(abs(vecs_smoothed(:,3:4))))]);