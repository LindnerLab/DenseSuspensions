function [Psmall, Plarge] = Verlet_compression_object(Psmall, Plarge, Rsmall, Rlarge, timestep, Vinit, Channel)
    Nsmall = size(Psmall,1);
    Nlarge = size(Plarge,1);
    Nparticles = Nsmall + Nlarge;
    dL_step = Channel(1)*0.5*0.01;
    Channel(1) = Channel(1)*1.5;
    
    for i = 1:Nparticles
        phi = rand()*2*pi - pi;
        if i <= Nsmall
            P(i) = ParticleVerlet([Psmall(i,1); Psmall(i,2)], Rsmall, [cos(phi);sin(phi)]*Vinit);
        else
            P(i) = ParticleVerlet([Plarge(i-Nsmall,1); Plarge(i-Nsmall,2)], Rlarge, [cos(phi);sin(phi)]*Vinit);
        end
    end
    
    for i = 1:10000
        for j = 1:Nparticles
            UpdatePosition(P(j),timestep);
            Pos = P;
            Pos(j) = [];
            UpdateVelocity(P(j), Channel, Pos);
        end
        
        if mod(i,100) == 0
            Channel(1) = Channel(1) - dL_step;            
        end
    end
    
    Psmall = [P(1:Nsmall).position]';
    Plarge = [P(Nsmall+1:end).position]';
end