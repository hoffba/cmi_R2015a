% ImageClass function
function stat = filtNOVA(self,varargin)
stat = false;
if self.check
    if isempty(license('inuse','Distrib_Computing_Toolbox'))
        TF = license('checkout','Distrib_Computing_Toolbox');
    else
        TF = true;
    end
    
    go = false;
    if nargin==1
        % Input filter parameters:
        prompt = {'Images to filter:','Radius of moving window cube',...
                  'Filter passes','Initial noise estimate',...
                  'Trapezoidal plateau (cannot be larger than max extent)',...
                  'Trapezoidal max extent'};
        title = 'NOVA filter parameters';
        default = {num2str(1:self.dims(4)),'3', '2', '200', '0', '3'};
        answer = inputdlg(prompt,title,1,default);
        if ~isempty(answer) && ~any(isnan(str2double(answer(2:end))))
            vec = str2num(answer{1});
            answer = str2double(answer(2:end));
            n = answer(1);
            run = answer(2);
            sDev0 = answer(3);
            p = answer(4);
            d = answer(5);
            go = true;
        end
    else
    end
    
    % Filter the selected images
    if go
        if TF
            timg = PnFiltNOVA(self.mat(:,:,:,vec),n,run,sDev0,p,d);
        else
            timg = nFiltNOVA(self.mat(:,:,:,vec),n,run,sDev0,p,d);
        end
        if ~isempty(timg)
            self.mat(:,:,:,vec) = timg;
            stat = true;
        end
    end
end
