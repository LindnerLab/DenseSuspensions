close all;
clear all;

directory = 'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\';
files = dir(fullfile(directory,'Preprocessed\V1\MAT','*.mat'));
nFiles = length(files);

for i = 1:nFiles
    Pall_struct = load(fullfile(files(i).folder,files(i).name)).file;
    Pall_CSV = [[Pall_struct.x]', [Pall_struct.y]', [Pall_struct.r]',ones(length(Pall_struct),1)*(i-1)];
    
    writematrix(Pall_CSV, fullfile(directory,'Preprocessed\V1\CSV',[sprintf('%05d',i),'.csv']),...
        'Delimiter', 'comma');
end