fn = '2016_06_07_0000.abf';
cd 'C:\Users\Administrator\Documents\Molecular Devices\pCLAMP\Data\'
[data,sampling_interval_us,header] = abfload([fn]);
Fs = 1/(sampling_interval_us*1e-6);

