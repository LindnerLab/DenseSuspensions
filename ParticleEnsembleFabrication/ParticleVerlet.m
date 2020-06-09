classdef ParticleVerlet < handle
    properties
        position
        dposition
        r
        velocity
    end
    
    methods
        function obj = ParticleVerlet(position, r, velocity)
            obj.position = position;
            obj.dposition = [0; 0];
            obj.r = r;
            obj.velocity = velocity;
        end
        function obj = UpdatePosition(obj, timestep)
            [obj.position] = [obj.position] + [obj.velocity]*timestep + [obj.dposition];
            [obj.dposition] = [0;0];
        end
        function obj = UpdateVelocity(obj, Channel, P)
            for i = 1:2
                if obj.position(i) < obj.r
                    obj.velocity(i) = -obj.velocity(i);
                    obj.dposition(i) = obj.r - obj.position(i);
                elseif obj.position(i) > Channel(i) - obj.r
                    obj.velocity(i) = -obj.velocity(i);
                    obj.dposition(i) = Channel(i) - obj.r - obj.position(i);
                end
            end
            
            Nparticles = length(P);
            d = sqrt(sum([[P.position]-[obj.position]].^2));
            for i = 1:Nparticles
                if d(i) < obj.r + P(i).r
                    %https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
                    d_vec = obj.velocity;
                    n_vec = [P(i).position-obj.position];
                    n_len = d(i);
                    
                    r_vec = d_vec - ((2*dot(d_vec,n_vec)/(n_len.^2))*n_vec);
                    obj.velocity = r_vec;
                    
                    obj.dposition = -0.5*(n_vec./n_len)*(obj.r+P(i).r-n_len);
                end
            end
        end
    end
end