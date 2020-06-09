%all units used are in mm, as the stage requires input in mm (and for the
%calculations it doesnt matter).

% close all;
% clear all;
verbose = true;

ChannelW = 5; %Width of the microfluidic channel to be filledin mm
ChannelL = 7.5; %Length of the microfluidic channel to be filled in mm
Phi_RCP = 0.66;
CR = 2/3;

Rsmall = 0.05; %Radius of the small particles in mm
Rlarge = 0.075; %Radius of the large particles in mm
AR = 1; %Ratio of Asmall over Alarge

timestep = 0.002;
Vinit = 1;
Nrecursion = 0;

[Psmall, Plarge] = Random_placement(ChannelW, ChannelL, Phi_RCP, Rsmall, Rlarge, AR, CR, Nrecursion);

if verbose
    Nsmall = length(Psmall);
    Nlarge = length(Plarge);
    
    figure;
    hold on;
    viscircles([Psmall(:,1), Psmall(:,2)],ones(Nsmall,1)*Rsmall);
    viscircles([Plarge(:,1), Plarge(:,2)],ones(Nlarge,1)*Rlarge);
    xlim([0 ChannelL]);
    ylim([0 ChannelW]);
end

[Psmall, Plarge] = Verlet_compression_vector(Psmall, Plarge, Rsmall, Rlarge, timestep, Vinit, [ChannelL, ChannelW],CR);

Psmall = ShortestPath(Psmall);
Plarge = ShortestPath(Plarge);


if verbose
    figure;
    viscircles([Psmall(:,1), Psmall(:,2)],ones(Nsmall,1)*Rsmall);
    viscircles([Plarge(:,1), Plarge(:,2)],ones(Nlarge,1)*Rlarge);
    xlim([0 ChannelL]);
    ylim([0 ChannelW]);
end

