function [mask] = CircularMask(r, Img_size)
    xCenter = r;                                    % X-position of the center of the object to look for. It is centered at (r,r), meaning that the peaks in the convoluted image will be at (x+r,y+r), instead of (x,y)!
    yCenter = r;                                    % Y-position of the center of the object to look for.
    xCircle = r*cos(0:pi/15:2*pi) + xCenter;        % X-coordinates of a circle centered at (xCenter, yCenter) with radius r.
    yCircle = r*sin(0:pi/15:2*pi) + yCenter;        % Y-coordinates of a circle centered at (xCenter, yCenter) with radius r.
    
    [X, Y] = meshgrid(1:Img_size(2),1:Img_size(1)); % Create a mesh of X and Y coordinates needed to create the mask (mesh should be equal in size to the image!)
    mask = inpolygon(X,Y,xCircle,yCircle);          % Check for each x and y in the mesh if it is inside or outside the 'object' on the mask.
end