clear all; close all; clc

% INPUT STIMULUS AND DATA FILES
fn1 = '2016_08_09_0004.abf';
fn2 = 'freqstim_f0.5to100_amp10_10sec_5kHz_N16.abf';
% fn2 = 'steps-ampn10to10-1sec-5kHz-N41.abf';
% fn2 = 'whitenoise_std1_Fs5kHz_30sec.abf';
% fn2 = 'freqstim_f0.5to100_amp50_5sec_5kHz_N27.abf';
% fn2 =  'freqstim_f0.5to100_amp25_10sec_5kHz_N16.abf';
% fn2 = 'steps-ampn10to10-1sec-5kHz-N41.abf';

% Import data
cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\'
[data,sampling_interval_us,header] = abfload([fn1]);
Fs = 1/(sampling_interval_us*1e-6);
dt = 1/Fs;

% Import stimulus waveform
cd 'C:\Users\Administrator\Desktop\Outfiles\'
[stimulus,sampling_interval_us_stim,header_stim] = abfload([fn2]);

N = input('N = ');

if N == 16
    freqstim = [0.5,1,2,4,6,8,10,15,20,25,30,35,40,45,50,100];
    disp(['Selected ' num2str(N) ' frequency values. Parsing...']);
elseif N == 20
    freqstim = [0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40];
    disp(['Selected ' num2str(N) ' frequency values. Parsing...']);
elseif N == 27
    freqstim = [0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100];
    disp(['Selected ' num2str(N) ' frequency values. Parsing...']);
else
    freqstim = input(['Freq (N=' num2str(N) '): ']);
    disp(['Selected ' num2str(N) ' frequency values. Parsing...']);
end

nonraw = 1:N;

% Define Ir, Xc 
% Because wave
stimyn = input('Stimulus? (y/n OR 1/0): ');

if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
if length(data) - length(stimulus) > 0
    data00 = data;data0=zeros(1,length(data));
    m=1;
    while length(data0) - length(stimulus) > 0 
        clear data0
        if m == 1
            indcut = round((length(data) - length(stimulus))/2);
        else
            indcut = round((length(data00) - length(stimulus)));
        end
        
        for j = 1:size(data,2)
            for k = 1:size(data,3)
                if m == 1
                    data0(:,j,k) = data(indcut:end-indcut,j,k);
                else
                    data0(:,j,k) = data00(1:end-indcut,j,k);
                end
            end
        end
        m=m+1;
    end
elseif length(stimulus) - length(data) > 0
    while length(stimulus) - length(data) > 0
        indcut = round((length(stimulus) - length(data)));
        stimulus = stimulus(1:end-indcut);
        data0 = data;
    end
else
    data0 = data;
end
else
    data0 = data;
end
% Partition the data
for k = 1:size(data,3)
for j = 1:N
    if N>1 && size(data,3) ==1
        xL = length(data0)/N;
        Ir0{j} = data0(1+(j-1)*xL:j*xL,k); Ir{j} = Ir0{j} - mean(Ir0{j});
        if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
            Xc0{j} = stimulus(1+(j-1)*xL:j*xL); Xc{j} = Xc0{j} - mean(Xc0{j});
        else
            Xc{j} = zeros(1,length(Ir{j}));
        end
        tvec{j} = 0:dt:(length(Ir{j})-1)*dt;
    elseif N == 1 && size(data,3) > 1
        xL = length(data0)/N;
        Ir0{k} = data0(:,k)'; Ir{k} = Ir0{k} - mean(Ir0{k});
        if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
            Xc0{k} = stimulus(1+(j-1)*xL:j*xL); Xc{k} = Xc0{k} - mean(Xc0{k});
        else
            Xc{k} = zeros(1,length(Ir{k}));
        end
        tvec{k} = 0:dt:(length(Ir{k})-1)*dt;
        nonraw = 1:size(data,3);
    elseif N > 1 && size(data,3) > 1
        xL = length(data0)/N;
        Ir0{j+(k-1)*N} = data0(1+(j-1)*xL:j*xL,k); Ir{j+(k-1)*N} = Ir0{j+(k-1)*N} - mean(Ir0{j+(k-1)*N});
        if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
            Xc0{j+(k-1)*N} = stimulus(1+(j-1)*xL:j*xL); Xc{j+(k-1)*N} = Xc0{j+(k-1)*N} - mean(Xc0{j+(k-1)*N});
        else
            Xc{j+(k-1)*N} = zeros(1,length(Ir{j+(k-1)*N}));
        end
        tvec{j+(k-1)*N} = 0:dt:(length(Ir{j+(k-1)*N})-1)*dt;
        nonraw = 1:N*size(data,3);
    end
end
end

disp('Saving...')
cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\'
mkdir(fn1(1:end-4)); cd(fn1(1:end-4));
save('Extracted Data.mat');
disp('Finished.')


%% Unsupervised Analysis
clear all; close all; clc;

date = '2016_08_09';        % INPUT
nanalyze = 6;               % INPUT

for j = 1:nanalyze
    if j-1 < 10
        fn = [date '_000' num2str(j-1)];
        disp(['Analyzing ' fn '...']);
        cd(['C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\' fn])
        spikesortinganalysis(pwd,0,0,1,0);
    elseif j-1 < 100
        fn = [date '_00' num2str(j-1)];
        disp(['Analyzing ' fn '...']);
        cd(['C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\' fn])
        spikesortinganalysis(pwd,0,0,1,0);
    end
end
disp('COMPLETE');

%% Reformat for wave_clus
clear all; close all; clc

% INPUT STIMULUS AND DATA FILES
cd 'C:\Users\Administrator\Documents\MATLAB\Simulator'
load('C_Easy2_noise02.mat');

fn1 = '2016_06_13_0002.abf';
cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\'
cd(fn1(1:end-4));
load('Extracted Data.mat');

samplingInterval = 1/Fs;
data = data';
chan = 1;
startData = 0;
par.sr = Fs;sr = Fs;
par.fmin_detect = 0.5*Fs;
par.fmax_detect = 0.1*Fs;
fmin_detect = par.fmin_detect;
fmax_detect = par.fmax_detect;


save('Extracted Data-wave_clus.mat');
disp('Finished.')

Get_spikes('Extracted Data-wave_clus.mat');
Do_clustering('Extracted Data-wave_clus_spikes.mat')
figfile = importdata('fig2print_Extracted Data-wave_clus.png');
close all
imshow(figfile)





