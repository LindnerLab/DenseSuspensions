function [Psmall, Plarge] = Particle_compression(Psmall, Plarge, Rsmall, Rlarge, timestep, Vg, ChannelL, ChannelW, pool)
    Nsmall = size(Psmall,1);
    Nlarge = size(Plarge,1);
    Nparticles = Nsmall + Nlarge;

    for i = 1:Nparticles
        if i <= Nsmall
            P(i) = Particle(Psmall(i,1), Psmall(i,2), Rsmall, Vg);
        else
            P(i) = Particle(Plarge(i-Nsmall,1), Plarge(i-Nsmall,2), Rlarge, Vg);
        end
    end

    count = 0;

    while count < 10000
        for i = 1:Nparticles
            UpdatePosition(P(i), timestep);
            IsOnBoundary(P(i), ChannelW);
            UpdatedV(P(i), P, i);
        end
        count = count + 1;
    end
    
    for i = 1:Nparticles
        Psmall = [[P(1:Nsmall).x]',[P(1:Nsmall).y]'];
        Plarge = [[P(Nsmall+1:end).x]',[P(Nsmall+1:end).y]'];
    end
end
