classdef Particle < handle
    properties
        x           %x position of the particle
        y           %y position of the particle
        r           %radius of the particle
        Vg          %gravitational force on particle
        dVx
        dVy
        Vx
        Vy
        XminBoundary = false;
        YminBoundary = false;
        YmaxBoundary = false;
    end
    
    methods
        function obj = Particle(x, y, r, Vg)
            obj.x = x;
            obj.y = y;
            obj.r = r;
            obj.Vg = Vg;
            obj.Vx = 0;
            obj.Vy = 0;
        end
        function obj = UpdatePosition(obj, timestep)
            obj.x = obj.Vx*timestep + obj.x;
            obj.y = obj.Vy*timestep + obj.y;
        end
        function obj = UpdatedV(obj, Positions, n)
            V = [-obj.Vg, 0];
            if obj.XminBoundary
                V = [0, 0];
            else
                Positions(n) = [];
                Nparticles = length(Positions);
                dx = [Positions.x] - obj.x;
                dy = [Positions.y] - obj.y;
                d = sqrt(dx.^2+dy.^2);
                for i = 1:Nparticles
                    if d(i) <= obj.r + Positions(i).r
                        theta = angle(dx(i) + 1i*dy(i));
                        phi = acos(dot([V(1), V(2)],[dx(i), dy(i)])/(sqrt(V(1)^2+V(2)^2)*sqrt(dx(i)^2+dy(i)^2)));
                        len = sin(phi)*sqrt(V(1)^2 + V(2)^2);
                        V(1) = cos(-theta)*len*sin(phi);
                        V(2) = sin(-theta)*len*sin(phi);
                        
                        obj.x = Positions(i).x - dx(i)*((obj.r + Positions(i).r)/d(i));
                        obj.y = Positions(i).y - dy(i)*((obj.r + Positions(i).r)/d(i));
                    end
                end
            end
            obj.Vx = V(1);
            obj.Vy = V(2);
        end
        function obj = IsOnBoundary(obj, ChannelW)
            if obj.x <= obj.r
                obj.XminBoundary = true;
                obj.x = obj.r;
            end
            
            if obj.y <= obj.r
                obj.YminBoundary = true;
                obj.y = obj.r;
            end
            
            if obj.y >= ChannelW - 2*obj.r
                obj.YmaxBoundary = true;
                obj.y = ChannelW - 2*obj.r;
            end
        end
    end
end