% ImageClass function
% Smooth image using Gaussian filter
function imgSmooth(self,tvec)
if self.check && (nargin==2) && all(tvec>0) && all(tvec<=self.dims(4))
    str = 3; % strength of filter
    prompt = {'Smooth dimension: (0=all,1=row,2=col,3=slc)','Filter strength:'};
    answer = inputdlg(prompt,'Smooth',1,{'0',num2str(str)});
    if ~isempty(answer)
        fdim = str2num(answer{1});
        if any((fdim<0) | (fdim>3))
            fdim = []; % invalid dimension index
        end
        str = str2double(answer{2});
        if ~isempty(fdim) && ~isnan(str) && (str>0)
            for iv = 1:length(tvec)
                if fdim % 1D filter
                    for i = 1:length(fdim)
                        % shift vdim to first dimension for filtering, then back
                        tmat = shiftdim(self.mat(:,:,:,tvec(iv)),vdim-1);
                        tmat = imfilter(tmat,fspecial('gaussian',[max(3,str),1],str));
                        self.mat(:,:,:,tvec(iv)) = shiftdim(tmat,mod(4-vdim,3));
                    end
                else % 3D filter
                    self.mat(:,:,:,tvec(iv)) = filtGaussSep(self.mat(:,:,:,tvec(iv)),str);
                end
            end
        end
    end
end