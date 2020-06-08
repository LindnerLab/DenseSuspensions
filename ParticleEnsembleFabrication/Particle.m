classdef Particle < handle
    properties
        x           %x position of the particle
        y           %y position of the particle
        r           %radius of the particle
        Fg          %gravitational force on particle
        Vx
        Vy
        ax
        ay
        Fx          %force in x-direction
        Fy          %force in y-direction
        XminBoundary = false;
        YminBoundary = false;
        YmaxBoundary = false;
    end
    
    methods
        function obj = Particle(x, y, r, Fg)
            obj.x = x;
            obj.y = y;
            obj.r = r;
            obj.Fg = Fg;
            obj.Vx = 0;
            obj.Vy = 0;
            obj.ax = 0;
            obj.ay = 0;
            obj.Fx = 0;
            obj.Fy = 0;
        end
        function obj = UpdatePosition(obj, timestep)
            obj.x = 0.5*obj.ax*timestep^2 + obj.Vx*timestep + obj.x;
            obj.y = 0.5*obj.ay*timestep^2 + obj.Vy*timestep + obj.y;
        end
        function obj = UpdateForce(obj, Positions, n)
            obj.Fx = obj.Fg;
            Positions(n) = [];
            Nparticles = length(Positions);
            dx = obj.x - Positions.x;
            dy = obj.y - Positions.y;
            d = sqrt(dx.^2-dy.^2);
            for i = 1:Nparticles
                if d <= obj.r + Positions(i).r
                    theta = angle(dx + 1i*dy);
                    
                end
            end
            
        end
        function obj = UpdateAcceleration(obj)
            obj.ax = obj.ax + obj.Fx;
            obj.ay = obj.ay + obj.Fy;
        end
        function obj = UpdateVelocity(obj, timestep)
            obj.Vx = obj.Vx + obj.ax*timestep;
            obj.Vy = obj.Vy + obj.ay*timestep;
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