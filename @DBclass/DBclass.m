classdef DBclass < handle
    % Image database class, for finding/opening known data
    %   Connects to SQLite database: clindata.db
    
    properties (SetAccess=private, GetAccess=public)
        tbCheck = false;    % Check whether 'Database_Toolbox' is available
        db_location = '';   % .db file location
        conn                % Object containing current database connection
    end
    
    methods
        % Contstructor w/ option to load specific database
        function self = DBclass(dbstr)
            % Input: dbstr - string containing location of database to load
            if license('test','Database_Toolbox')
                self.tbCheck = true;
                if nargin
                    if ~(ischar(dbstr) && exist(dbstr,'file'))
                        warning('Invalid database input.')
                    else
                        self.db_location = dbstr;
                    end
                else
                    % Set to default location:
                    netdrive = findNetDrives('maize');
                    self.db_location = fullfile(netdrive,'Database','clindata.db');
                end
                if ~isempty(self.db_location)
                    tobj = sqlite(self.db_location);
                    if isconnection(tobj)
                        self.conn = tobj;
                    end
                end
            end
        end
        function stat = connect(self,dbstr)
            stat = false;
            if nargin==1
                dbstr = self.db_location;
            elseif ~(ischar(dbstr) && exist(dbstr,'file') && strcmp(dbstr(end-2:end),'.db'))
                warning('Invalid input.');
                return;
            end
            self.conn = sqlite(dbstr);
        end
        function addMetric(self,name,descr,units)
            if ischar(name) && ischar(descr) && ischar(units)
                insert(self.conn,'metrics',{'name','description','units'},{name,dscr,units});
            else
                error('Invalid inputs.');
            end
        end
        function addMeasurement(self,metric,imageID,val)
            metnames = fetch(self.conn,sprintf('SELECT name FROM metrics WHERE name = %s'));
            imnums = fetch(self.conn,'SELECT imageid FROM images');
            if ischar(metricID) && ismember(metric,metnames) ...
                    && isnumeric(imageID) && isscalar(imageID) && ismember(imageID,imnums) ...
                    && isnumeric(val) && isscalar(val)
                insert(self.conn,'metrics',{'name','description','units'},{name,dscr,units});
            else
                error('Invalid inputs.');
            end
        end
        function addImage(self,tag,
            
        end
        function addScan
        end
        function addSubject
        end
        function addStudy
        end
        function subj2study
        end
        function catalogDICOM(self)
        end
        % Add image to database
        function addImg(self)
        end
        function delete(self)
        end
    end
    
end

