close all;
clearvars -except Img;

directory = 'F:\Lars\Oscillatory Compression\20200713 Soft Particle\Avg100_Amp80_Per120_Back25\';
files = dir(fullfile(directory,'Preprocessed','*.mat'));
nFiles = length(files);

verbose = true;
Save_files = true;

maxD = 10;
param = struct('mem',10,...
    'dim',2,...
    'good',1000,...
    'quiet',0);

%% Below this point no user input needed

positions = [];
for i = 1:nFiles
    load(fullfile(files(i).folder,files(i).name));
    nParticles = length(file);
    positions = [positions;[file(:).x]',[file(:).y]',[file(:).r]',ones(nParticles,1)*i];
end

positions_tracked = track(positions,maxD,param);

if verbose
    figure
    imshow(Img);
    hold on;
    posx = NaN(2000,1384);
    posy = NaN(2000,1384);
    post = NaN(2000,1384);
    
    for i = 1:max(positions_tracked(:,5))
        idx = positions_tracked(:,5) == i;
        posx(i,1:sum(idx)) = positions_tracked(idx,1);
        posy(i,1:sum(idx)) = positions_tracked(idx,2);
        post(i,1:sum(idx)) = positions_tracked(idx,4);
%         plot(positions_tracked(idx,1),positions_tracked(idx,2));
    end
   
    posx_norm = posx - (max(posx,[],2,'omitnan')+min(posx,[],2,'omitnan'))./2;
    posx_norm = posx_norm./max(posx_norm,[],2,'omitnan');
    
    posy_norm = posy - (max(posy,[],2,'omitnan')+min(posy,[],2,'omitnan'))./2;
    posy_norm = posy_norm./max(posy_norm,[],2,'omitnan');
    
    figure;
    hold on;
    for i = 101:200
        plot(post(i,:),posx_norm(i,:));
    end
    
    figure;
    hold on;
    for i = 100:102
        plot(post(i,:),posy_norm(i,:));
    end
    
    xpos = positions_tracked(idx,1);
    ypos = positions_tracked(idx,2);
    t = positions_tracked(idx,4);
    
    xpos = xpos - mean(xpos);
    ypos = ypos - mean(ypos);
    
    xpos = xpos./max(xpos);
    ypos = ypos./max(ypos);
    
    figure
    hold on;
    plot(t,xpos);
    plot(t,ypos);
    xlabel('t (s)');
    ylabel('(pos - mean) / Amplitude');
    legend('x','y');
    
    figure
    histogram(mod(positions(:,1),1),100);
    xlabel('remainder x-pos mod(1) (px)');
    ylabel('# occurrences');
    
    for i = 3
        n = sum(not(isnan(posx(i,:))));
        p = plot(posx(i,:),posy(i,:));
        cd = [uint8(jet(n)*255) uint8(ones(n,1))].';

        drawnow
        set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
        xlabel('x-pos (px)');
        ylabel('y-pos (px)');
    end
    
    
end




for i = 1:nFiles
    load(fullfile(files(i).folder,files(i).name));
    pos = [[file(:).x]',[file(:).y]',[file(:).r]'];
    writematrix(pos,fullfile(directory,'Preprocessed','CSV',[sprintf('%05d',i),'.csv']));    
end








