function m = curve8(type)
    im=zeros(5,5);
    im(1,1) = 1;
    im(2,1) = 1;
    im(3,1) = 1;
    im(4,1) = 1;
    im(5,1) = 1;
    [x,y] = contour2d(im,type);
    m = [y,x];         
end