classdef RegClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % GUI and options
        guicheck = true; % Determines whether to show GUI or not
        h                % structure of all GUI handles
        
        % Elastix parameters / functions:
        elxObj  % ElxClass object for storing/executing Elastix commands
        ind     % Selected schedule step
        odir = '';    % Elastix/Transformix Output directory
        
        % Warping Extras:
        jac     = false; % Spatial |Jacobian| map
        jacmat  = false; % Jacobian Matrices
        def     = false; % Deformation Maps
        
        % Image interface
        cmiObj           % set of two CMIclass objects (constructor)
                         %       (1) Reference Image
                         %       (2) Homologous Image
        points = cell(1,2); % {1 x 2}[n x 3] control points in spatial coordinates
        hpts = nan(1,2); % handles to point plots
        T0defVal = 0;      % Default image value for the OneShot
        
        % Pre-processing Options:
        clamp = [ -inf , inf ; -inf , inf ];
        ftype = 'Gaussian';
        filtN = zeros(1,3);
        dilateN = zeros(1,3);
        
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
        function self = RegClass
            self.cmiObj = CMIclass;
            self.cmiObj(2) = CMIclass;
            self.elxObj = ElxClass;
            
            % Default queue file location:
            self.qfile = fullfile(fileparts(which('cmi')),'elastix','qfile.txt');
            
            %main Callback set function
            self.initMain;
        end
        % RegClass destructor
        function delete(self)
            delete(self.cmiObj(isvalid(self.cmiObj)));
            delete(self.elxObj);
            if isfield(self.h,'regFig') && ishandle(self.h.regFig)
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