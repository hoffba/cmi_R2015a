function nMinf = OptimAffine( refImg, homImg, pts )
%function to be converted into an anonymous function for optimization
%process
refpts = pts(:,:,1);
hompts = pts(:,:,2);
nImg = afftform(homImg, 1, 'linear', refpts, hompts);
nImg(isnan(nImg)) = min(nImg(:));
nMinf = -mutualinfo(nImg, refImg);

%display mutual information and histogram per every iteration
disp(strcat('\n current MI is ', num2str(-nMinf)));

end

