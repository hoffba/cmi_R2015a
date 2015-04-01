function G = sobel3D(img)

d = size(img);
if length(d)==3
    G = zeros(size(img));
    h = [ 1,2,1 ; 2,4,2 ; 1,2,1 ];
    h = cat(3,h,zeros(3),-h);
    for i = 1:3
        G = sqrt(G.^2 + convn(img,h,'same').^2);
        h = shiftdim(h,1);
    end
end