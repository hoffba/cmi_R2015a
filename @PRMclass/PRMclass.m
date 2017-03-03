classdef PRMclass < Mat3Dclass
    properties (SetAccess=private, GetAccess=public)
        % Inherited properties
        %   mat         % Matrix of indexed PRM values
        %   dims        % PRM matrix dimensions
        %   check       % Check that matrix is available
        
        % Current PRM options
        thresh = zeros(0,4);   % [n x 4] matrix : row = [vec1 vec2 m b]
                               %       evaluated: d2 > m*d1 + b
        cutoff = zeros(0,3);   % [n x 3] matrix : row = [vec min max]
        
        % Pre-filtering for PRM:
        %   ** Default to [3x3] Median filter
        filtchk = true;
        filttype = 'median';
        filtstr = [3,3];
        
        cmap     % RGB columns
        prmmap   % {n x 2} cell array:
                 %       col1 = [n x nthresh] logical array for matching to C
                 %       col2 = label of region
                 %       row # --> PRM region index
        nprm     % Number of PRM colors
        dlabels  % cell array of image labels
        dvec     % vector of image #'s in use, corresponds to dlabels
                 %      dvec(1) is the current/dynamic image #
                 %      dvec(2) is the baseline image
                 %      dvec(3:end) are other requested images
                 
        % Scatterplot options
        SPopts = struct('Xvec',{nan},'Yvec',{nan},'Xmin',{nan},'Ymin',{nan},...
                        'Xmax',{nan},'Ymax',{nan},'show',{true},'Nmax',{5000});
        mask     % MaskClass associated with PRM selections
                 
        prmdir   % Directory storing PRM default settings
        
        guicheck = true;       % Flag to display results
        hfscatter               % Handle to PRM scatterplot figure
        hascatter               % Handle to PRM scatterplot axes
        hsscatter               % Handle to PRM scatterplot object
        htscatter               % Handle to PRM scatterplot stats textbox
        
    end
    methods
        % Constructor
        function self = PRMclass(obj)
            if ~nargin
                obj = [];
            end
            self = self@Mat3Dclass(obj);
            % Find files containing PRM settings
            self.prmdir = fullfile(fileparts(which('cmi')),'PRMdefs');
            if ~isdir(self.prmdir)
                mkdir(self.prmdir);
            end
            pdefs = dir(fullfile(self.prmdir,'*PRMdef*.mat'));
            pdefs = {pdefs(:).name};
            idef = find(~cellfun(@isempty,strfind(pdefs,'RGB')),1);
            if isempty(idef)
                idef = 1;
            end
            self.loadPRMdefs(fullfile(self.prmdir,pdefs{idef}));
            self.nprm = size(self.cmap,1);
                tval = self.thresh(:,1:2);
                tval = unique(tval(:));
            self.dvec = tval(:)';
            self.dlabels = cellstr(strcat('Dim',num2str((1:length(self.dvec))')));
        end
        % Deleter method
        function delete(self)
            if ~isempty(self.hfscatter) && ishandle(self.hfscatter)
                close(self.hfscatter);
            end
        end
    end
end