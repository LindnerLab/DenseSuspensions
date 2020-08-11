function [mask] = CircInCircMask(Img_size,r)
    mask_temp = zeros(Img_size(1),Img_size(2),3);
    xCenter = r(1);                                             % X-position of the center of the object to look for. It is centered at (r,r), meaning that the peaks in the convoluted image will be at (x+r,y+r), instead of (x,y)!
    yCenter = r(1);                                             % Y-position of the center of the object to look for.
    for i = 1:3
        xCircle = r(i)*cos(0:pi/15:2*pi) + xCenter;             % X-coordinates of a circle centered at (xCenter, yCenter) with radius r.
        yCircle = r(i)*sin(0:pi/15:2*pi) + yCenter;             % Y-coordinates of a circle centered at (xCenter, yCenter) with radius r.
        [X, Y] = meshgrid(1:Img_size(2),1:Img_size(1));
        
        mask_temp(:,:,i) = inpolygon(X,Y,xCircle,yCircle);      % Check for each x and y in the mesh if it is inside or outside the 'object' on the mask.
    end
    
    mask = mask_temp(:,:,1) - 5*mask_temp(:,:,2) + 5*mask_temp(:,:,3);
end