close all;
clear all;

directory = 'F:\Lars\Oscillatory Compression\20200630 Ensemble Oscillation\Ensemble_oscillation_compression\Avg130_Amp100_Per120\';
files = dir(fullfile(directory, 'RawData','*.tif'));
nFiles = length(files);
verbose = true;
Rsmall = [31, 28];

rect = [121,163,2183,2026];

% [~,rect] = imcrop(imread(fullfile(files(1).folder,files(1).name)));

SelectionCriteria = struct('Property',[],'Value',[],'Criteria',[]);
% The particle shouldn't be too small
SelectionCriteria(1).Property = 'Area';
SelectionCriteria(1).Value = pi*(0.8*Rsmall(1))^2;
SelectionCriteria(1).Criteria = 'Greater';

% The particle shouldn't be too large either
SelectionCriteria(2).Property = 'Area';
SelectionCriteria(2).Value = pi*(1.2*Rsmall(1))^2;
SelectionCriteria(2).Criteria = 'Smaller';

% The particle should be roughly round (0 is perfectly round, 1 is a line)
SelectionCriteria(3).Property = 'Eccentricity';
SelectionCriteria(3).Value = 0.4;
SelectionCriteria(3).Criteria = 'Smaller';

for i = 1%:nFiles
    ImgRaw = double(imcrop(imread(fullfile(files(i).folder,files(i).name)),rect));
    ImgSize = size(ImgRaw);
    ImgBg = ImgRaw - 2E4;
    ImgBg(ImgBg < 0) = 0;
    ImgGauss = imgaussfilt(ImgBg,ImgSize(1)/10);
    ImgGauss(ImgGauss < 2000) = 2000;
    ImgCorr = ImgBg-ImgGauss;
    ImgCorr(ImgBg == 0) = 0;
    
    mask = AnnulusMask(Rsmall, ImgSize);
    [Psmall] = FindParticlesConvolution(ImgCorr,Rsmall(1),SelectionCriteria,0.55,mask,verbose);
end
