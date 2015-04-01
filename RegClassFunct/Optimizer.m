function nMinf = Optimizer(refImg, homImg, F)
%Basic transform and cost function for creating a new image using mutual
%information cost function
B = F(:,:,1);
t = F(:,:,2);
T = B*t;
homImg = im3affine(homImg, T, 1);
homImg(isnan(homImg)) = min(homImg(:));
nMinf = -mutualinfo(homImg, refImg);
disp(strcat('current mutual information is ', num2str(-nMinf)));
end


