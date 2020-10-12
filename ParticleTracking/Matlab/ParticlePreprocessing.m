close all;
clear all;

tic;

directory = 'F:\Lars\Oscillatory Compression\20200820 Soft Particle Oscillation\Avg75_Amp50_Per120_Back25\';
version = 'V1';
files = dir(fullfile(directory,'RawData',version,'_00001.tif'));
output = strcat('Preprocessed\',version);
nFiles = length(files);

load_settings = true;
verbose = true;
Inv_img = true; 
Dev_gauss = true;
Save_files = false;  
Parallel_processing = false;

if load_settings
    Settings = open(fullfile(directory,output,'Settings.mat')).Settings;
    Rsmall = Settings.Rsmall;
    Rlarge = Settings.Rlarge;
    Inv_img = Settings.Inv_img;
    thresh = Settings.thresh;
    Crop = Settings.Crop;
    SelectionCriteria = Settings.SelectionCriteria;
else 
    Rsmall = [19, 22];
    Rlarge = [26, 30];

    [Img,Crop] = imcrop(imread(fullfile(files(1).folder,files(1).name)));
    
    thresh = 3.2E4;
    
    SelectionCriteria = struct('Property',[],'Value',[],'Criteria',[]);

    % The convoluted peak corresponding to a particle shouldn't be too small
    % (size depends on the size of the counted region of the mask)
    SelectionCriteria(1).Property = 'Area';
    SelectionCriteria(1).Value = 35;
    SelectionCriteria(1).Criteria = 'Greater';

    % The convoluted peak shouldn't be too large either
    % (size depends on the size of the counted region of the mask)
    SelectionCriteria(2).Property = 'Area';
    SelectionCriteria(2).Value = 180;
    SelectionCriteria(2).Criteria = 'Smaller';
end

%% Find Particles -- Below this no input required
if ~exist(fullfile(directory,'Preprocessed'),'dir')
    mkdir(fullfile(directory,'Preprocessed'));
    mkdir(fullfile(directory,'Preprocessed',version));
    mkdir(fullfile(directory,'Preprocessed',version,'MAT'));
    mkdir(fullfile(directory,'Preprocessed',version,'CSV'));
    
elseif ~exist(fullfile(directory,'Preprocessed',version),'dir')
    mkdir(fullfile(directory,'Preprocessed',version));
    mkdir(fullfile(directory,'Preprocessed',version,'MAT'));
    mkdir(fullfile(directory,'Preprocessed',version,'CSV'));
end

if Parallel_processing && isempty(gcp('nocreate'))
        p = parpool(10);
end

for i = 1:nFiles
    if not(isfile(fullfile(directory,output,[sprintf('%05d',i),'.mat']))) || not(Save_files)
        Img = imcrop(imread(fullfile(files(i).folder,files(i).name)),Crop);
        Img_size = size(Img);
        
        if Dev_gauss
            Img = double(Img)./imgaussfilt(double(Img),50);
            Img = uint16(Img.*(39050/1.0245));
        end
        
        if Inv_img
            Img_inv = imcomplement(Img);
            Img_thresh = Img_inv;
        else
            Img_thresh = Img;
        end

        Img_thresh(Img_thresh < thresh) = 0;

        mask = AnnulusMask(Rsmall, Img_size);
        [Psmall] = FindParticlesConvolution(Img_thresh,mean(Rsmall),SelectionCriteria,0.2,mask,verbose);
        Nsmall = length(Psmall);

        mask = AnnulusMask(Rlarge, Img_size);
        [Plarge] = FindParticlesConvolution(Img_thresh,mean(Rlarge),SelectionCriteria,0.2,mask,verbose);
        Nlarge = length(Plarge);

        %% Remove false positives
        Psmall = Psmall(Psmall(:,1) > Rsmall(1),:);
        Plarge = Plarge(Plarge(:,1) > Rlarge(1),:);

        Nsmall = length(Psmall);
        Nlarge = length(Plarge);

        Distances = zeros(Nlarge,Nsmall);
        for j = 1:Nsmall
            Distances(:,j) = sqrt((Plarge(:,1)-Psmall(j,1)).^2+(Plarge(:,2)-Psmall(j,2)).^2);
        end

        Distances_min = min(Distances);

        Psmall = Psmall(Distances_min > Rsmall(1)+Rlarge(1),:);
        Nsmall = length(Psmall);

        if verbose
            figure
            imshow(Img);
            viscircles([Psmall(:,1), Psmall(:,2)],ones(Nsmall,1)*mean(Rsmall),...
                'EdgeColor','b',...
                'LineWidth',1.5);
            viscircles([Plarge(:,1), Plarge(:,2)],ones(Nlarge,1)*mean(Rlarge),...
                'EdgeColor','r',...
                'LineWidth',1.5);
        end

        %% Save found particle locations
        if Save_files
            Pall = struct('x',[],'y',[],'r',[]);
            X = [Psmall(:,1);Plarge(:,1)];
            Y = [Psmall(:,2);Plarge(:,2)];
            R = [ones(Nsmall,1)*mean(Rsmall);ones(Nlarge,1)*mean(Rlarge)];

            for j = 1:Nsmall+Nlarge
                Pall(j).x = X(j);
                Pall(j).y = Y(j);
                Pall(j).r = R(j);
            end
            SaveParallel(fullfile(directory,output,'MAT',[sprintf('%05d',i),'.mat']),...
                         fullfile(directory,output,'CSV',[sprintf('%05d',i),'.csv']),...
                         Pall);
        end
    end
end

if Save_files
    Settings = struct('Rsmall',Rsmall,...
                      'Rlarge',Rlarge,...
                      'Inv_img',Inv_img,...
                      'thresh',thresh,...
                      'Crop',Crop,...
                      'SelectionCriteria',SelectionCriteria);
    save(fullfile(directory,output,'Settings.mat'),'Settings');
end

toc





















