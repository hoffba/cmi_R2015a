% Calculates Mutual Information of two images of the same dimensions/size
function I = cmi_mi(v1,v2,nbin)
I = [];
if (nargin==2) && (length(v1)==length(v2))
    % Get range of values
    r1 = [min(v1),max(v1)];
    r2 = [min(v2),max(v2)];
    % Generate 1D histograms
    [n1, bin1] = histc(x,xedges);
    [n2, bin2] = histc(y,yedges);

    %xbin, ybin zero for out of range values 
    % (see the help of histc) force this event to the 
    % first bins
    xbin(xbin == 0) = inf;
    ybin(ybin == 0) = inf;

    if xnbin >= ynbin
        xy = ybin*(xnbin) + xbin;
          indexshift =  xnbin; 
    else
        xy = xbin*(ynbin) + ybin;
          indexshift =  ynbin; 
    end

    xyuni = unique(xy);
    xyuni(end) = []; 
    hstres = histc(xy,xyuni);
    clear xy;

    histmat = zeros(ynbin,xnbin);
    histmat(xyuni-indexshift) = hstres;

    
end