function imgout = cmi_filter(img,n,varargin)
% cmi_filter(img,k,'Name',Value,...)
% Image filter for 1D, 2D, or 3D
% Inputs: img = image matrix for filtering
%         n   = filter neighborhood (radius)

p = inputParser;
addRequired(p,'Img',@(x)isnumeric(x)||islogical(x));
addRequired(p,'N',@isnumeric);
addParamValue(p,'Shape','circ');
addParamValue(p,'Type','gauss',@ischar);
parse(p,img,n,varargin{:});
opts = p.Results;

% Kernel properties
nN = length(opts.N);
N = zeros(1,3);
N(1:nN) = opts.N;

switch opts.Type
    case {'median','erode','dilate','open','close'}
        switch opts.Type
            case 'erode'
                tType = {'min'};
            case 'dilate'
                tType = {'max'};
            case 'open'
                tType = {'max','min'};
            case 'close'
                tType = {'min','max'};
            otherwise
                tType = {opts.Type};
        end
        imgout = img;
        for i = 1:length(tType)
            imgout = ordfilt3(imgout,N,'Type',tType{i});
        end
    case {'average','gauss','unsharp'}
        h = ones(2*N+1);
        if any(strcmp(opts.Type,{'gauss','unsharp'}))
                sig = (2*N+1)/2/2.354;
                [x,y,z] = ndgrid(-N(1):N(1),-N(2):N(2),-N(3):N(3));
                h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        end
        h = h.*cmi_kmask(N,opts.Shape);
        h = h/sum(h(:));
        imgout = imfilter(img,h);
        if strcmp(opts.Type,'unsharp')
            imgout = img - imgout*0.8;
        end
    case 'wiener' % only implemented for 2D
        imgout = img;
        for i = 1:size(img,3)
            imgout(:,:,i) = wiener2(img(:,:,i),N(1:2));
        end
end




