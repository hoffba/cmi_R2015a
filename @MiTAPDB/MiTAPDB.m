classdef MiTAPDB < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        dbsourcename = 'GalbanPDB';
        dbname = 'galban_prod';
        serv = 'galbanpdb.med.umich.edu';
        portno = 3330;
        conn
    end
    methods
        function delete(self)
            if isa(self.conn,'database.jdbc.connection') && isopen(self.conn)
                close(self.conn);
            end
        end
        function stat = connect(self)
            if isa(self.conn,'database.jdbc.connection') && isopen(self.conn)
                stat = true;
            else
                % List all available DB sources
                list = listDataSources();
                if ismember(self.dbsourcename,list.Name)
                    % DB connection is already set up
                    opts = databaseConnectionOptions(self.dbsourcename);
                else
                    % Need to set up DB connection 
                    % First find Connector/J installation
                    [fname,fpath] = uigetfile('*.jar');
                    if ~fname
                        fprintf(['Download the Connector/J release at: https://dev.mysql.com/downloads/connector/j/\n',...
                                 '    (Typicall installation location: C:\\Program Files (x86)\\MySQL\\\n)']);
                        return;
                    else
                        opts = databaseConnectionOptions("jdbc","MySQL");
                        opts = setoptions(opts,...
                            'DataSourceName',self.dbsourcename,...
                            'JDBCDriverLocation',fullfile(fpath,fname),...
                            'DatabaseName','galban_prod',...
                            'Server',self.serv,...
                            'PortNumber',self.portno,...
                            'useSSL','true',...
                            'serverTimezone','America/New_York');
                        saveAsDataSource(opts);
                    end
                end
                
                % User input for Lvl-2
                [uname,pw] = uilogin('Title','Lvl-2 Login');
                
                % Test Connection
                stat = testConnection(opts,uname,pw);
    
                if stat
                    self.conn = database(self.dbsourcename,uname,pw);
                else
                    fprintf('Connection test failed. Removing from database source list.\n');
                    deleteDataSource(self.dbsourcename);
                end
            end
        end
        function ind = addStudy(self,name,descr)
            if self.connect
                % Make sure there isn't a row with the same name
                data = select(self.conn,sprintf('SELECT idStudy FROM Studies WHERE Name LIKE "%s"',name));
                if isempty(data)
                    T = table();
                    T.Name = {name};
                    T.Description = {descr};
                    sqlwrite(self.conn,'Studies',T);
                    ind = select(self.conn,'SELECT COUNT(*) FROM Studies');
                    ind = ind.COUNT___;
                else
                    fprintf('Study "%s" already exists.\n',name);
                    ind = data.idStudy;
                end
            end
        end
        function ind = addSubject(self,name)
            if self.connect
                ind = 0;
                if length(name)>10
                    warning('Subject name is too long. Maximum 10 characters.');
                else
                    data = select(self.conn,['SELECT Name FROM Subjects WHERE Name LIKE ',name]);
                    if isempty(data)
                        T = table();
                        T.Name = {name};
                        sqlwrite(self.conn,'Subjects',T);
                        ind = true;
                    else
    
                        fprintf('Subject %s already exists.\n',name);
                    end
                end
            end
        end
        function stat = addDrive(self,name,WinStr,UnixStr)
            if self.connect
                stat = false;
                T = sqlread(self.conn,'Drives');
                if ismember(name,T.Name)
                    fprintf('Drive %s already exists.\n',name);
                else
                    T = table();
                    T.Name = {name};
                    T.WinStr = {WinStr};
                    T.UnixStr = {UnixStr};
                    sqlwrite(self.conn,'Drives',T);
                    stat = true;
                end
            end
        end
        function stat = addMetric(self,name,descr,units)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addMeasurement(self,imgID,segID,segVal,metricStr,val)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addImage(self,tagID,driveID,scanID,procID,loc)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addProcess(self,name,descr)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addTag(self,name,descr,modality)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addScan(self)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function stat = addSegmentation(self)
            stat = false;
            T = sqlread(self.conn,'Metrics');
        end
        function n = tableCount(self,Tname)
            if self.connect
                T = select(self.conn,['SELECT COUNT(*) FROM ',Tname]);
                n = T.COUNT___;
            end
        end
        function T = selectFromTable(self,Tname,varargin)
            if self.connect

                qry = ['SELECT * FROM ',Tname,'WHERE' ];

                T = select(self.conn,['SELECT * FROM ',Tname,'WHERE']);
            end
        end
    end
end

