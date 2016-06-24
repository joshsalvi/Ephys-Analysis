clear all; close all; clc

% INPUT STIMULUS AND DATA FILES
fn1 = '2016_06_22_0003.abf';
fn2 = 'P841-freqstim_f0.1to100log_amp10_30sec_1kHz_N25.abf';

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
stimyn = input('Stimulus? (y/n): ');

if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
if length(data) - length(stimulus) > 0
    indcut = round((length(data) - length(stimulus))/2);
    data0 = data(indcut:end-indcut);
else
    data0 = data;
end
else
    data0 = data;
end
% Partition the data
for j = 1:N
    xL = length(data0)/N;
    Ir0{j} = data0(1+(j-1)*xL:j*xL); Ir{j} = Ir0{j} - mean(Ir0{j});
    if stimyn == 1 || stimyn == 'Y' || stimyn == 'y'
        Xc0{j} = stimulus(1+(j-1)*xL:j*xL); Xc{j} = Xc0{j} - mean(Xc0{j});
    else
        Xc{j} = zeros(1,length(Ir{j}));
    end
    tvec{j} = 0:dt:(length(Ir{j})-1)*dt;
end

disp('Saving...')
cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\'
mkdir(fn1(1:end-4)); cd(fn1(1:end-4));
save('Extracted Data.mat');
disp('Finished.')


%% Unsupervised Analysis
clear all; close all; clc;

date = '2016_06_22';        % INPUT
nanalyze = 6;               % INPUT

for j = 1:nanalyze
    if j-1 < 10
        fn = [date '_000' num2str(j-1)];
        disp(['Analyzing ' fn '...']);
        cd(['C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\' fn])
        spikesortinganalysis(pwd,0,1,1,1);
    elseif j-1 < 100
        fn = [date '_00' num2str(j-1)];
        disp(['Analyzing ' fn '...']);
        cd(['C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\MATLAB\' fn])
        spikesortinganalysis(pwd,0,1,1,1);
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





