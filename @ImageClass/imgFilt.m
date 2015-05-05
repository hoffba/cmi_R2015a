% ImageClass function
% Performs 2D image filter
function imgFilt(self,tvec,pars)
if self.check && (nargin>=2) && (tvec>0) && (tvec<=self.dims(4))
    go = false;
    if nargin<3
        % Ask for user inputs:
        answer = inputdlg({'Image Vector(s)','Filter Neighborhood',...
                            'Filter Type ([w]iener,[m]edian,[g]aussian,[u]nsharp mask)'},...
                          'Median Filter',1,...
                          {num2str(tvec),'3','m'});
        if ~isempty(answer)
            tvec = str2num(answer{1});
            fstr = str2num(answer{2});
            ftype = strsplit(answer{3},' ');
            go = true;
        end
    elseif iscell(pars) && (length(pars)==3)
        tvec = pars{1};
        fstr = pars{2};
        ftype = pars{3}; % cell array of strings
        go = true;
    end
    if go
        nv = length(tvec);
        nstr = length(fstr);
        nt = length(ftype);
        if (nstr==1) && (nv>1)
            fstr = fstr*ones(1,nv);
        end
        if (nt==1) && (nv>1)
            ftype = repmat(ftype,1,nv);
        end
        for iv = 1:nv
            timg = self.mat(:,:,:,tvec(iv));
            ifilt = fstr(iv)*[1,1];
            switch ftype{iv}
                case 'm'
                    func = @(x)medfilt2(x,ifilt);
                case 'w'
                    func = @(x)wiener2(x,ifilt);
                case 'g'
                    ifilt = round(fstr(iv)*6);
                    func = @(x)imfilter(x,fspecial('gaussian',ifilt,fstr(iv)));
                case 'u'
                    func = @(x)imsharpen(x,'Radius',fstr(iv)*6);
                otherwise
                    func = [];
            end
            if ~isempty(func)
                hw = waitbar(0,'Applying 2D filter ...');
                for islc = 1:self.dims(3)
                    timg(:,:,islc) = feval(func,timg(:,:,islc));
                    waitbar(islc/self.dims(3),hw);
                end
                self.mat(:,:,:,tvec(iv)) = timg;
                delete(hw);
            end
        end
    end
end