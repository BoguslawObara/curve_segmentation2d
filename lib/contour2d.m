function [x,y] = contour2d(im,type) 
%%  contour2d - find contour points of the curve
%
%   INPUT:
%       im      - binary image with contour/curve.
%       type    - type of the contour: 'close' or 'open'.
%
%   OUTPUT:
%       x,y     - x and y coordinaties of the contour
%
%   AUTHOR: 
%       Boguslaw Obara
        
[row, col] = find(im,1);
contour = bwtraceboundary(im, [row, col], 'N');
if strcmp(type,'open')
    contour = contour(1:round(size(contour,1)/2),:);
elseif strcmp(type,'close')
    contour = contour;
end
y = contour(:,1); 
x = contour(:,2); 

end