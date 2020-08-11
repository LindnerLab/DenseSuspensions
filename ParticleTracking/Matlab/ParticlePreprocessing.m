close all;
clear all;

directory = 'F:\Lars\Oscillatory Compression\20200713 Soft Particle\Avg100_Amp80_Per120_Back25\';
files = dir(fullfile(directory,'RawData','*.tif'));
nFiles = length(files);

load_settings = true;
verbose = false;
Inv_img = true;
Save_files = true;  
Parallel_processing = true;

if load_settings
    uiopen('*.mat');
    Rsmall = Settings.Rsmall;
    Rlarge = Settings.Rlarge;
    Inv_img = Settings.Inv_img;
    thresh = Settings.thresh;
    Crop = Settings.Crop;
    SelectionCriteria = Settings.SelectionCriteria;
else 
    Rsmall = [19, 22];
    Rlarge = [26, 30];

    [Img,Crop] = imcrop(imread(fullfile(directory,'RawData',files(1).name)));
    
    thresh = 3.6E4;
    
    SelectionCriteria = struct('Property',[],'Value',[],'Criteria',[]);

    % The convoluted peak corresponding to a particle shouldn't be too small
    % (size depends on the size of the counted region of the mask)
    SelectionCriteria(1).Property = 'Area';
    SelectionCriteria(1).Value = 50;
    SelectionCriteria(1).Criteria = 'Greater';

    % The convoluted peak shouldn't be too large either
    % (size depends on the size of the counted region of the mask)
    SelectionCriteria(2).Property = 'Area';
    SelectionCriteria(2).Value = 200;
    SelectionCriteria(2).Criteria = 'Smaller';
end

%% Find Particles -- Below this no input required
if Parallel_processing && isempty(gcp('nocreate'))
        p = parpool(6);
end

parfor i = 1:nFiles
    if not(isfile([directory,'Preprocessed\',sprintf('%05d',i),'.mat'])) || not(Save_files)
        Img = imcrop(imread(fullfile(directory,'RawData',files(i).name)),Crop);
        Img_size = size(Img);

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
            viscircles([Psmall(:,1), Psmall(:,2)],ones(Nsmall,1)*mean(Rsmall),'EdgeColor','b');
            viscircles([Plarge(:,1), Plarge(:,2)],ones(Nlarge,1)*mean(Rlarge),'EdgeColor','r');
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
            SaveParallel([directory,'Preprocessed\',sprintf('%05d',i),'.mat'],Pall);
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
    save(fullfile(directory,'Preprocessed','Settings.mat'),'Settings');
end





















