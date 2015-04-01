classdef DBclass < handle
    % Image database class, for finding/opening known data
    %   Detailed explanation goes here
    
    properties (SetAccess=private, GetAccess=public)
        tbCheck = false;    % Check whether 'Database_Toolbox' is available
        dbConnObj           % Object containing current database connection
        user = 'tech';      % Database username
        upw = 'tech';       % Database password
    end
    
    methods
        % Contstructor w/ option to load specific database
        function self = DBclass(dbstr)
            % Input: dbstr - string containing location of database to load
            if license('test','Database_Toolbox')
                self.tbCheck = true;
                if ~(nargin && ischar(dbstr) && exist(dbstr,'file'))
                    % Ask user to select database
                    [fname,fpath] = uigetfile('*.mdb','Find Database');
                    if fname
                        dbstr = fullfile(fpath,fname);
                    else
                        dbstr = [];
                    end
                else
                    dbstr = [];
                end
                if ~isempty(dbstr)
                    tobj = database(dbstr,self.user,self.upw);
                    if isconnection(tobj)
                        self.dbConnObj = tobj;
                    end
                end
            end
        end
        % Add image to database
        function addImg(self)
        end
    end
    
end

