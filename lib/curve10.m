function m = curve10(type)
    im=zeros(10,10);
    im(1,1) = 1;
    im(2,1) = 1;
    im(3,2) = 1;
    im(4,2) = 1;
    im(5,2) = 1;
    im(6,2) = 1;
    im(7,3) = 1;
    im(8,3) = 1;
    im(9,3) = 1;
    im(10,4) = 1;
    im(11,4) = 1;
    im(12,4) = 1;
    im(13,5) = 1;        
    [x,y] = contour2d(im,type);
    m = [y,x];         
end