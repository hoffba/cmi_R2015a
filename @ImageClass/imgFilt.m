% ImageClass function
% Performs 2D image filter
function imgFilt(self,vec,ftype,opts)
% Inputs:
%   vec = image number to filter
%   ftype = (string) type of filter to use
%   opts = (vector) filter parameters

filts = {'Median','Gauss','Wiener','Average','Unsharp','AniDiffWavelet','ThreshWavelet','Mode'};
if self.check && (nargin>=2) && all(vec>0) && all(vec<=self.dims(4))
    go = true;
    
    % Determine filter type:
    if nargin<3
        ftype = listdlg('ListString',filts,'SelectionMode','single',...
            'Name','Filter Type:','PromptString','What type of filter?');
        if ~isempty(ftype)
            ftype = filts{ftype};
            go = true;
        else
            return;
        end
    elseif ~ischar(ftype) || ~any(strcmpi(ftype,filts))
        error('ImageClass::imgFilt - Invalid filter type input.')
    end
    
    % Determine filter options:
    switch lower(ftype)
        case 'median'
            str = {'Neighborhood'};
            defs = {'3 3'};
            nchk = true;
            func = { @(x,opts) medfilt2(x,opts{1}) ,...
                     @(x,opts) medfilt3(x,opts{1}) };
        case 'gauss'
            str = {'Sigma'};
            defs = {'0.5 0.5'};
            nchk = true;
            func = { @(x,opts) imgaussfilt(x,opts{1}) ,...
                     @(x,opts) imgaussfilt3(x,opts{1}) };
        case 'wiener'
            str = {'Neighborhood'};
            defs = {'3 3'};
            nchk = true;
            func = @(x,opts) wiener2(x,opts{1});
        case 'average'
            str = {'Neighborhood'};
            defs = {'3 3'};
            nchk = true;
            func = @(x,opts) imfilter(x,ones(opts{1})/prod(opts{1}));
        case 'unsharp'
            str = {'Sigma','Amount','Threshold'};
            defs = {'1','0.8','0.7'};
            nchk = true(1,3);
            func = @(x,opts) imsharpen(x,'Radius',opts{1},'Amount',opts{2},'Threshold',opts{3});
        case 'anidiffwavelet'
            str = {'Wavelet Filter (wfilters)','Wavelet Levels',...
                '(Diff) Iterations','(Diff) dt','(Diff) Conduction',...
                'Conduction Function (1/2)'};
            defs = {'sym8','1','50','1/7','50','1'};
            nchk = [false,true(1,5)];
            func = @(x,opts) anidiffWT(x,opts{1},opts{2},opts{3},opts{4},opts{5},opts{6});
        case 'threshwavelet'
            str = {'Wavelet Filter (wfilters)','Wavelet Levels','Threshold','Hard/Soft (0/1)'};
            defs = {'sym8','1','20','1'};
            nchk = [false,true(1,3)];
            func = @(x,opts) threshWT(x,opts{1},opts{2},opts{3},opts{4});
        case 'mode'
            str = {'Neighborhood'};
            defs = {'3 3'};
            nchk = true;
            func = @(x,opts) colfilt(x,opts{1},'sliding',@mode);
    end
    
    % User input for filter options:
    if (nargin<4)
        answer = inputdlg([{'Image(s):'},str],ftype,1,[{num2str(vec)},defs]);
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
        hw = waitbar(0,'Applying filter:');
        tmat = nan([self.dims(1:3),length(vec)]);
        for v = 1:length(vec)
            if iscell(func) && numel(opts{1})==3
                waitbar(0,hw,'Applying 3D image filter:');
                tmat(:,:,:,v) = feval(func{2},self.mat(:,:,:,vec(v)),opts);
%                 tmat(:,:,:,v) = medfilt3(self.mat(:,:,:,vec(v)),opts{1});
            else
                for i = 1:self.dims(3)
                    tmat(:,:,i,v) = feval(func{1},self.mat(:,:,i,vec(v)),opts);
                    ct = ct+1;
                    waitbar(ct/ntot,hw,['Applying 2D image filter: ',num2str(ct)]);
                end          
            end
        end
        self.mat(:,:,:,vec) = tmat;
        delete(hw);
    end
    
end

