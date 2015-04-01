function H = histFunct(refImg, homImg)
%HISTFUNCT function that outputs histogram given two images along MI info

xvals = refImg;
yvals = homImg;
binlenX = range(xvals)/(1+log2(numel(xvals)));
binlenY = range(yvals)/(1+log2(numel(yvals)));
binVecX = 1:binlenX:range(xvals);
binVecY = 1:binlenY:range(yvals);
H = hist2(xvals,yvals,binVecX, binVecY);
end

