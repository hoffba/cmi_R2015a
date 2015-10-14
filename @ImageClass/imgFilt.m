% ImageClass function
% Performs 2D image filter
function imgFilt(self,vec,ftype,opts)
% Inputs:
%   vec = image number to filter
%   ftype = (string) type of filter to use
%   opts = (vector) filter parameters

filts = {'Median','Gauss','Wiener','Average','Unsharp','AniDiffWT','Mode'};
if self.check && (nargin>=2) && all(vec>0) && all(vec<=self.dims(4))
    go = false;
    
    % Determine filter type:
    if nargin<3
        ftype = listdlg('ListString',filts,'SelectionMode','single',...
            'Name','Filter Type:','PromptString','What type of filter?');
        if ~isempty(ftype)
            ftype = filts{ftype};
            go = true;
        end
    elseif ~ischar(ftype) || ~any(strcmpi(ftype,filts))
        error('ImageClass::imgFilt - Invalid filter type input.')
    end
    
    % User input for filter options:
    if (nargin<4)
        N = [];
        switch lower(ftype)
            case 'median'
                opts = {'Neighborhood'};
                defs = {'3 3'};
                nchk = true;
                func = @(x,opts) medfilt2(x,opts{1});
            case 'gauss'
                opts = {'Sigma'};
                defs = {'0.5 0.5'};
                nchk = true;
                func = @(x,opts) imgaussfilt(x,opts{1});
            case 'wiener'
                opts = {'Neighborhood'};
                defs = {'3 3'};
                nchk = true;
                func = @(x,opts) wiener2(x,opts{1});
                N = zeros(1,self.dims(3));
            case 'average'
                opts = {'Neighborhood'};
                defs = {'3 3'};
                nchk = true;
                func = @(x,opts) imfilter(x,ones(opts{1}));
            case 'unsharp'
                opts = {'Sigma','Amount','Threshold'};
                defs = {'1','0.8','0.7'};
                nchk = true(1,3);
                func = @(x,opts) imsharpen(x,'Radius',opts{1},'Amount',opts{2},'Threshold',opts{3});
            case 'anidiffwt'
                opts = {'Wavelet Filter (wfilters)','Levels','Iterations','Conduction'};
                defs = {'dmey','4','20','70'};
                nchk = [false,true(1,3)];
                func = @(x,opts) anidiffWT(x,opts{1},opts{2},opts{3},1/7,opts{4},1);
            case 'mode'
                opts = {'Neighborhood'};
                defs = {'3 3'};
                nchk = true;
                func = @(x,opts) colfilt(x,opts{1},'sliding',@mode);
        end
        answer = inputdlg([{'Image(s):'},opts],ftype,1,[{num2str(vec)},defs]);
        if isempty(answer)
            go = false;
        else
            vec = str2num(answer{1});
            opts = answer(2:end);
            opts(nchk) = cellfun(@(x)sscanf(x,'%f')',answer([false,nchk]),'UniformOutput',false);
            go = true;
        end
    elseif ~iscell(opts)
        error('ImageClass::imgFilt - Invalid filter options, must be cell array of vectors.')
    end
    
    % Perform the filtering:
    if go
        ct = 0;
        ntot = self.dims(3)*length(vec);
        hw = waitbar(0,'Applying 2D image filter:');
        for v = 1:length(vec)
            for i = 1:self.dims(3)
                val{:} = feval(func,self.mat(:,:,i,vec(v)),opts);
                self.mat(:,:,i,vec(v)) = val{1};
                if ~isempty(N)
                    N(i) = val{2};
                end
                ct = ct+1;
                waitbar(ct/ntot,hw,['Applying 2D image filter: ',num2str(ct)]);
            end
        end
        delete(hw);
        
        if ~isempty(N)
            figure,plot(1:self.dims(3),N);
            title('Wiener Filter Noise Estimates');
            xlabel('Slice Number');
            ylabel('Image Noise');
        end
    end
    
end



% if self.check && (nargin>=2) && all(tvec>0) && all(tvec<=self.dims(4))
%     go = false;
%     if nargin<3
%         % Ask for user inputs:
%         answer = inputdlg({'Image Vector(s)','Filter Neighborhood',...
%             'Filter Type ([a]verage,[w]iener,[m]edian,[g]aussian,[u]nsharp mask,anisotropic wavelet [d]iffusion)'},...
%             'Filter Image:',1,...
%             {num2str(tvec),'3','m'});
%         if ~isempty(answer)
%             tvec = str2num(answer{1});
%             fstr = str2num(answer{2});
%             ftype = strsplit(answer{3},' ');
%             go = true;
%         end
%     end
%     if go && isnumeric(fstr) && (iscellstr(ftype) || ischar(ftype))
%         if ischar(ftype)
%             ftype = {ftype};
%         end
%         nv = length(tvec);
%         nstr = length(fstr);
%         nt = length(ftype);
%         if (nstr==1) && (nv>1)
%             fstr = fstr*ones(1,nv);
%         end
%         if (nt==1) && (nv>1)
%             ftype = repmat(ftype,1,nv);
%         end
%         for iv = 1:nv
%             timg = self.mat(:,:,:,tvec(iv));
%             ifilt = fstr(iv)*[1,1];
%             switch ftype{iv}
%                 case 'a'
%                     func = @(x)imfilter(x,ones(fstr));
%                 case 'm'
%                     func = @(x)medfilt2(x,ifilt);
%                 case 'w'
%                     func = @(x)wiener2(x,ifilt);
%                 case 'g'
%                     ifilt = round(fstr(iv)*6);
%                     func = @(x)imfilter(x,fspecial('gaussian',ifilt,fstr(iv)));
%                 case 'u'
%                     func = @(x)imsharpen(x,'Radius',fstr(iv)*6);
%                 otherwise
%                     func = [];
%             end
%             if ~isempty(func)
%                 hw = waitbar(0,'Applying 2D filter ...');
%                 for islc = 1:self.dims(3)
%                     timg(:,:,islc) = feval(func,timg(:,:,islc));
%                     waitbar(islc/self.dims(3),hw);
%                 end
%                 self.mat(:,:,:,tvec(iv)) = timg;
%                 delete(hw);
%             end
%         end
%     end
% end