classdef RegClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % GUI and options
        guicheck = true; % Determines whether to show GUI or not
        h                % structure of all GUI handles
        
        % Elastix parameters / functions:
        elxObj  % ElxClass object for storing/executing Elastix commands
        ind     % Selected schedule step
        odir = '';    % Elastix/Transformix Output directory
        waitchk = false; % default to run independently
        
        % Warping Extras:
        jac     = false; % Spatial |Jacobian| map
        jacmat  = false; % Jacobian Matrices
        def     = false; % Deformation Maps
        
        % Transform Extras:
        Tvoi    = false; % Transform moving VOI
        Tsurf   = false; % Transform fixed VOI surface mesh points
        
        % Image interface
        cmiObj           % set of two CMIclass objects (constructor)
                         %       (1) Reference Image
                         %       (2) Homologous Image
        points = cell(1,2); % {1 x 2}[n x 3] control points in spatial coordinates
        hpts = nan(1,2); % handles to point plots
        T0defVal = 0;      % Default image value for the OneShot
        
        % Pre-processing Options:
        clamp = [ -inf , inf ; -inf , inf ];
        ftype = 'Median';
        filtN = zeros(2,3);
        dilateN = zeros(2,3);
        unmaskval = nan(1,2);
        histeq = false;
        
        % Transformix options:
        Tfx = struct('out',{''},'par',{''},'jac',{false},'jacmat',{false},'def',{false},...
            'nn',{false(0,1)},'fnames',{cell(0,1)});
        
        % Batch job object for dynamic queue
        qfile   % File containing queue of system commands for batch processing
        qnum = 1; % Number of batch jobs running the queue
        job     % Batch job object(s)
        
    end
    methods
        % RegClass constructor
        function self = RegClass(guicheck)
            
            % Default queue file location:
            self.qfile = fullfile(fileparts(which('cmi')),'elastix','qfile.txt');
           
            % Validate inputs:
            if nargin
                if ischar(guicheck)
                    switch guicheck
                        case 'qfile'
                            if exist(self.qfile,'file')
                                fid = fopen(self.qfile);
                                if fid>0
                                    go = true;
                                    ct = 0;
                                    while go
                                        tline = fgets(fid);
                                        if tline==-1
                                            fprintf('Number of registrations in queue = %u\n',ct);
                                            go = false;
                                        else
                                            disp(tline);
                                            ct = ct+1;
                                        end
                                    end
                                    fclose(fid);
                                else
                                    fprintf('Could not open qfile.txt\n');
                                end
                            else
                                fprintf('Qfile does not exist.\n');
                            end
                    end
                    delete(self);
                    return;
                elseif ~isempty(guicheck)
                    self.guicheck = logical(guicheck);
                end
            end
            self.cmiObj = CMIclass(self.guicheck);
            self.cmiObj(2) = CMIclass(self.guicheck);
            self.elxObj = ElxClass;
            
            
            %main Callback set function
            self.initMain;
        end
        % RegClass destructor
        function delete(self)
            if ~isempty(self.cmiObj) && isa(self.cmiObj,'CMIclass')
                delete(self.cmiObj(isvalid(self.cmiObj)));
            end
            if ~isempty(self.elxObj) && isa(self.elxObj,'ElxClass')
                delete(self.elxObj);
            end
            if isfield(self.h,'regFig') && ishghandle(self.h.regFig)
                delete(self.h.regFig)
            end
        end
        % GUI close function
        function GUIcloseFcn(self,~,~)
            answer = questdlg('Would you like to close the Reg program?','Close Reg');
            if strcmp(answer,'Yes')
                delete(self);
            end
        end
        % Set ButtonDownFcn on both CMIclass overlay images, used for placing points
        function setBDFcn(self,~,eventdata)
            if (nargin==3) && isa(eventdata,'event.PropertyEvent')
                set(eventdata.AffectedObject.hiover,'ButtonDownFcn',@self.addPoint);
            end
        end
        % Listener for updating points list
        function listenPoints(self,~,~)
            eopts = {'on','off'};
            str = eopts{1+isempty(self.points{1})};
            set(self.h.list_ref,'String',num2str(self.points{1},'% .2f '),...
                               'Enable',str,...
                               'Value',max(size(self.points{1},1),1));
            str = eopts{1+isempty(self.points{2})};
            set(self.h.list_hom,'String',num2str(self.points{2},'% .2f '),...
                               'Enable',str,...
                               'Value',max(size(self.points{2},1),1));
        end
        % Listener for change in Initial Transform
        function listenT0(self,~,~)
            if isfield(self.elxObj.Tx0,'Transform')
                str = self.elxObj.Tx0.Transform;
                val = 1;
            else
                str = '[None]';
                val = 0;
            end
            set(self.h.text_currTform,'String',str,'Value',val);
        end
    end
end