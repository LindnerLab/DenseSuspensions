function [Psmall, Plarge] = Verlet_compression_vectorV2(Psmall, Plarge, Rsmall, Rlarge, timestep, Vinit, Channel, CR)
    Nsmall = size(Psmall,1);
    Nlarge = size(Plarge,1);
    Nparticles = Nsmall + Nlarge;
    dL_step = 0.5*Rsmall;
    Nsteps = ceil((Channel(1)*(1-CR))/dL_step);
    Tstep = 1000;
    
    P = [Psmall;Plarge];
    R = [ones(Nsmall,1)*Rsmall;ones(Nlarge,1)*Rlarge];
    phi = rand(Nparticles,1)*2*pi - pi;
    V = [cos(phi),sin(phi)]*Vinit; 
    dP = zeros(Nparticles,2);
    
    for i = 1:Nsteps*Tstep + 10*Tstep
        P = P + V*timestep + dP;
        dP = zeros(Nparticles,2);
        if mod(i,100) == 1
            neighbors = NaN(Nparticles,30);
        end

        for j = 1:Nparticles
            Vtemp = V(j,:);
            dPtemp = zeros(2,1);
            for k = 1:2
                if P(j,k) < R(j)
                    Vtemp(k) = -Vtemp(:,k);
                    dPtemp(k) = R(j) - P(j,k);
                elseif P(j,k) > Channel(k) - R(j)
                    Vtemp(k) = -Vtemp(:,k);
                    dPtemp(k) = Channel(k) - R(j) - P(j,k);
                end
            end
            dP(j,:) = dPtemp;
            V(j,:) = Vtemp;
            
            if mod(i,100) == 1
                d = sqrt(sum((P-P(j,:)).^2,2));
                Nneighbor = 1;
                for k = 1:Nparticles
                    if d(k) < 6*Rsmall && not(k == j)
                        neighbors(j,Nneighbor) = k;
                        Nneighbor = Nneighbor + 1;
                    end
                end
            end
            
            Pother = P(neighbors(j,~isnan(neighbors(j,:))),1:2);
            Rother = R(neighbors(j,~isnan(neighbors(j,:))));
            d = sqrt(sum((Pother - P(j,:)).^2,2));
            for k = 1:sum(~isnan(neighbors(j,:)))
                if d(k) < R(j) + Rother(k)
                    %https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
                    d_vec = V(j,:);
                    n_vec = Pother(k,:)-P(j,:);
                    n_len = d(k);

                    r_vec = d_vec - ((2*dot(d_vec,n_vec)/(n_len.^2))*n_vec);
                    V(j,:) = r_vec;

                    dPtemp = -0.5*(n_vec./n_len)*(R(j)+Rother(k)-n_len);
                end
            end
            dP(j,:) = dPtemp;
        end
        
        if mod(i,Tstep) == 0 && i < Tstep*Nsteps
            Channel(1) = Channel(1) - dL_step;            
        end
    end
    
    Psmall = P(1:Nsmall,:);
    Plarge = P(Nsmall+1:end,:);
end