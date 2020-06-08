function [Psmall, Plarge] = Random_placement(ChannelW, ChannelL, Phi_RCP, Rsmall, Rlarge, AR)
Aparticles = ChannelW*ChannelL*Phi_RCP; %Determine area used by particles, given a Phi_RCP

Nsmall = ceil(Aparticles*sqrt(AR)/(2*pi*Rsmall^2));
Nlarge = ceil(Aparticles/(sqrt(AR)*2*pi*Rlarge^2));

Nparticles = Nsmall + Nlarge; %Total number of particles is the number of small particles plus the number of large particles.
Particles = [ones(Nsmall,1)*Rsmall; ones(Nlarge,1)*Rlarge]; %Makes a list of length Nparticles with Nsmall number of Rsmall entries and Nlarge number of Rlarge entries

Psmall = [];
Plarge = [];
Pparticles = [];
Dsmall = [];
Dlarge = [];

for i = Nparticles:-1:1
    found = false;
    Ntry = 1;
    idx = randi(i);
    r = Particles(idx);
    
    if i == Nparticles
        x = rand*ChannelL;
        y = rand*ChannelW;
        
        if r == Rsmall
            Psmall = [Psmall; [x,y]];
            Pparticles = [Pparticles; [x,y]];
        else
            Plarge = [Plarge; [x,y]];
            Pparticles = [Pparticles; [x,y]];
        end
        found = true;
    end
    
    while not(found) && Ntry < 1000
        x = rand*(ChannelL*1.5-2*r)+r;
        y = rand*(ChannelW-2*r)+r;
        
        if not(isempty(Psmall))
            Dsmall = [sqrt((Psmall(:,1)-x).^2+(Psmall(:,2)-y).^2)];
        end
        
        if not(isempty(Plarge))
            Dlarge = [sqrt((Plarge(:,1)-x).^2+(Plarge(:,2)-y).^2)];
        end
        
        if sum([Dsmall < Rsmall+r; Dlarge < Rlarge+r]) == 0
            Pparticles = [Pparticles; [x,y]];
            if r == Rsmall
                Psmall = [Psmall; [x,y]];
            else
                Plarge = [Plarge; [x,y]];
            end
            Particles = [Particles(1:idx-1,:); Particles(idx+1:end,:)];
            found = true;
        end
        Ntry = Ntry + 1;
    end
end

























