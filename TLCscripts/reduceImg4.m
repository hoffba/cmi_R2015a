function y = reduceImg4(img)
% average pixels 4x4
% 8/14/10 TC

[d1 d2] = size(img);
y = zeros(d1/4,d2/4);

y(:,:) = (img(1:4:d1,1:4:d2) + img(2:4:d1,1:4:d2) + img(3:4:d1,1:4:d2) + img(4:4:d1,1:4:d2) ...
        + img(1:4:d1,2:4:d2) + img(2:4:d1,2:4:d2) + img(3:4:d1,2:4:d2) + img(4:4:d1,2:4:d2) ...
        + img(1:4:d1,3:4:d2) + img(2:4:d1,3:4:d2) + img(3:4:d1,3:4:d2) + img(4:4:d1,3:4:d2) ...
        + img(1:4:d1,4:4:d2) + img(2:4:d1,4:4:d2) + img(3:4:d1,4:4:d2) + img(4:4:d1,4:4:d2) )/16;
           
clear img;