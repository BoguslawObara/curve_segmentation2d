function m = curveCircle2(type)
    s = 50;
    im = zeros(s,s);
    im(20,20) = 1;
    imd = imdilate(im,strel('disk',17));
    ime = imerode(imd,strel('disk',1));
    im = imsubtract(imd,ime);
    [x,y] = contour2d(im,type);
    m = [y,x];         
end