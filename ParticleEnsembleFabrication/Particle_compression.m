%function [Psmall, Plarge] = Particle_compression(Psmall, Plarge, Rsmall, Rlarge, timestep, Fg)
Psmall = [1,1];
Plarge = [];
timestep = 0.1;
Frep = 0.5;

Nsmall = length(Psmall);
Nlarge = length(Plarge);
Nparticles = Nsmall + Nlarge;

for i = 1:1
    if i <= Nsmall
        P(i) = Particle(Psmall(i,1), Psmall(i,2), Rsmall, -Fg);
    else
        P(i) = Particle(Plarge(i-Nsmall,1), Plarge(i-Nsmall,2), Rlarge, -Fg);
    end
end

cont = true;

while cont
    for i = 1
        UpdatePosition(P(i), timestep);
        IsOnBoundary(P(i), ChannelW);
        UpdateForce(P(i), P, i, Frep);
        UpdateAcceleration(P(i));
        UpdateVelocity(P(i), timestep);
        test = [P(i).x, P(i).y];
    end
end
