function y = reduceImg2(img)
% average pixels 2x2
% 8/14/07 HG

[d1 d2] = size(img);
y = zeros(d1/2,d2/2);

y(:,:) = (img(1:2:d1,1:2:d2) + img(2:2:d1,1:2:d2) ...
        + img(1:2:d1,2:2:d2) + img(2:2:d1,2:2:d2))/4;
           
clear img;