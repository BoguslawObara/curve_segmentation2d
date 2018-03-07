function m = curve9(type)
    im=zeros(5,5);
    im(1,1) = 1;
    im(2,2) = 1;
    im(3,3) = 1;
    im(4,4) = 1;
    im(5,5) = 1;
    [x,y] = contour2d(im,type);
    m = [y,x];         
end