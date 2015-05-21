% MaskClass function
% Perform morphological function on VOI
function morph(self,meth,options)
% meth: 'dilate' or 'erode'
% options: vector of morphological options
if (nargin == 3) && ~isempty(options) && isnumeric(options)
    opts = {'dilate','erode'};
    if any(strcmpi(meth,opts(1:2))) % Dilate / Erode
        % options: [xy-radius z-radius]
        %          * if z-rad = 0, makes 2D
        r = abs(options([1,1,2]));
        
        % Determine iterations (max single dilate/erode of 5)
        ni = max(ceil(r/5));
        
        hw = waitbar(0,['Performing ',meth,' (',num2str(r'),') - ',num2str(ni),' iterations']);
        for i = 1:ni
            rt = min(r,5);
            r = r - rt;
            se = bwellipsoid(rt);
            if strcmpi(meth,'dilate')
                self.mat = imdilate(self.mat,se);
            else
                self.mat = imerode(self.mat,se);
            end
            waitbar(i/ni,hw);
        end
        delete(hw);
    end
end