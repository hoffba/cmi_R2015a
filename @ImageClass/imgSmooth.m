% ImageClass function
% Smooth image using Gaussian filter
function imgSmooth(self,tvec)
if self.check && (nargin==2) && all(tvec>0) && all(tvec<=self.dims(4))
    str = 3; % strength of filter
    prompt = {'Smooth dimension: (0=all,1=row,2=col,3=slc)','Filter Type: 1=box, 2=Gaussian','Filter strength:'}; % added 'Filter Type...' 12Oct2015 CJG
    answer = inputdlg(prompt,'Smooth',1,{'0','1',num2str(str)}); % added '1'. 12Oct2015 CJG
    if ~isempty(answer)
        fdim = str2num(answer{1});
        if any((fdim<0) | (fdim>3))
            fdim = []; % invalid dimension index
        end
        str = str2double(answer{3}); % changed this from 2 to 3. 12Oct2015 CJG
        if answer{2}==2 % added 12Oct2015 CJG
            str_filt='gaussian';
        else
            str_filt='box';
        end
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
%                     self.mat(:,:,:,tvec(iv)) =
%                     filtGaussSep(self.mat(:,:,:,tvec(iv)),str); %
%                     commented out by CJG 12Oct2015
                    self.mat(:,:,:,tvec(iv)) = smooth3(self.mat(:,:,:,tvec(iv)),str_filt,[str str str]); % Added in by CJG 12Oct2015
                end
            end
        end
    end
end