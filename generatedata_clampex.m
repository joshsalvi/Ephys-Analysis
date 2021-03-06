function [tvec outvec] = generatedata_clampex(duration_sec,Fs,modtype,modparams,outfile)
%
% [tvec outvec] = generatedata_clampex(duration_sec,Fs,modtype,modparams,outfile)
%
% duration_sec: duration in seconds
% Fs: sampling rate
% modtype: 1,2,... defines type of modulation
% modparams: cell array with parameters for that modtype
% outfile: string containing output filename
%
% Coder:    Joshua D. Salvi
% Year:     2016
%

dt = 1/Fs;
duration_sec = duration_sec;
scalfac = [5.66 6.87 7.9 8.96 9.56 9.82];
cmdpzt = [0.1 0.5 1 2 3 4];
SF = 1/mean(scalfac);
    
if modtype == 1
    ampl = modparams{1};
    freq = modparams{2};
    for j = 1:length(freq)
        %tvec(:,j) = (0+(j-1)*(dt+Fs*dt*duration_sec)):dt:(Fs*dt*duration_sec+(j-1)*(dt+Fs*dt*duration_sec));
        tvec(:,j) = linspace(0+duration_sec*(j-1),duration_sec*j,Fs*duration_sec);
        outvec(:,j) = SF*ampl.*sin(2*pi*freq(j).*tvec(:,j));
    end
elseif modtype == 2
    ampl = modparams{1};
    freq = modparams{2};
    for j = 1:length(ampl)
        %tvec(:,j) = (0+(j-1)*(dt+Fs*dt*duration_sec)):dt:(Fs*dt*duration_sec+(j-1)*(dt+Fs*dt*duration_sec));
        tvec(:,j) = linspace(0+duration_sec*(j-1),duration_sec*j,Fs*duration_sec);
        outvec(:,j) = SF*ampl(j).*sin(2*pi*freq.*tvec(:,j));
    end
elseif modtype == 3
    ampl = modparams{1};
    for j = 1:length(ampl)
        %tvec(:,j) = (0+(j-1)*(dt+Fs*dt*duration_sec)):dt:(Fs*dt*duration_sec+(j-1)*(dt+Fs*dt*duration_sec));
        tvec(:,j) = linspace(0+duration_sec*(j-1),duration_sec*j,Fs*duration_sec);
        outvec(:,j) = SF*ampl(j).*randn(1,length(tvec));
    end
elseif modtype == 4
    ampl = modparams{1};
    delay = modparams{2};
    duration = modparams{3};
    for j = 1:length(ampl);        
        tvec(:,j) = linspace(0+duration_sec*(j-1),duration_sec*j,Fs*duration_sec);
        outvec(:,j) = zeros(1,length(tvec));
        qr1 = findnearest(tvec(:,1),delay);
        qr2 = findnearest(tvec(:,1),duration+delay);
        outvec(qr1:qr2,j) = ampl(j);
    end
end

sizeT = size(tvec);
tvec = reshape(tvec,sizeT(1)*sizeT(2),1).*1000;       % ms
outvec = reshape(outvec,sizeT(1)*sizeT(2),1);   % V

if exist('outfile')
    disp('Saving...')
    %m = [tvec,outvec];
    m = [outvec];
    dlmwrite(outfile,m,'\t')
    disp('Saved.')
end


end