
function job = runBatchFromQ(fname,varargin)
% Start / run batch Elastix queue from specified file:
% Inputs:   fname = File name of queue

job = [];
if nargin<2
    % This code starts the batch
    disp(['Starting batch for processing queue: ',fname])
    job = batch(@runBatchFromQ,1,{fname,1});
    disp(' ... done')
elseif varargin{1}==1
    
    go = true;
    while go
        % Grab next command from queue file:
        tname = fullfile(fileparts(fname),'temp.txt');
        ifid = fopen(fname,'rt');
        str = fgetl(ifid);

        % Re-save Q file without this line:
        ofid = fopen(tname,'wt');
        while ~feof(ifid)
           l = fgetl(ifid);
           fprintf(ofid,'%s\n',l);
        end
        fclose(ifid);
        fclose(ofid);
        copyfile(tname,fname,'f');
        delete(tname);

        if isempty(str) || (str(1)==-1)
            % If no more commands, delete the Q file and stop the loop
            delete(fname);
            go = false;
        else
            % Start process in new xterm:
            disp(['Processing system command: ',str]);
            system(str);
        end
    end
    
end
