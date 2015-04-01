classdef ElastixQueue < handle
    properties
        % File containing queue of Elastix calls
        % Set in constructor for ../cmi/elastix/qfile.txt
        qfile
        elxdir
    end
    methods
        function self = ElastixQueue
            self.qfile = fullfile(fileparts(which('cmi')),'elastix','qfile.txt');
            if ispc
            else
                self.elxdir = '/opt/elastix/bin';
            end
        end
        function enqueue(self,ffname,mfname,pfnames,outdir,varargin)
            % Syntax: @ElastixQueue.enqueue(ffname,mfname,pfnames,outdir,'Name',Value,...)
            % Required inputs (in order):
            %   ffname  = fixed image
            %   mfname  = moving image
            %   pfnames = Elastix parameter file (string or cell array of strings)
            %   outdir  = output directory
            % Optional inputs (Name/Value pairs):
            %   fMask   = fixed image mask
            %   mMask   = moving image mask
            %   t0      = initial transform parameter file
            %   threads = number of cores to use
            
            % Parse inputs:
            p = inputParser;
            isMHDfun = @(x) ischar(x) && (length(x)>4) && strcmp(x(end-3:end),'.mhd');
            isTXTfun = @(x) (ischar(x) && (length(x)>4) && strcmp(x(end-3:end),'.txt'));
            Pfun = @(x) (ischar(x) && (length(x)>4) && strcmp(x(end-3:end),'.txt')) ...
                         || (iscellstr(x) && all(cellfun(@length,x)>4) ...
                         && all(cellfun(@(xx)strcmp(xx(end-3:end),'.txt'),x)));
            addRequired(p,'f',isMHDfun);
            addRequired(p,'m',isMHDfun);
            addRequired(p,'p',Pfun);
            addRequired(p,'out',@(x) ischar(x) && isdir(x));
            addParamValue(p,'fMask','',isMHDfun);
            addParamValue(p,'mMask','',isMHDfun);
            addParamValue(p,'t0','',isTXTfun);
            addParamValue(p,'threads','',@(x) isinteger(x) && (length(x)==1) && (x>0));
            addParamValue(p,'orig','',isMHDfun);
            addParamValue(p,'jac',false,@(x) islogical(x) && (length(x)==1));
            addParamValue(p,'CleanUp',false,@(x) islogical(x) && (length(x)==1));
            parse(p,ffname,mfname,pfnames,outdir,varargin{:});
            
            % Generate Elastix command:
            pfnames = p.Results.p;
            if ischar(pfnames)
                pfnames = {pfnames};
            end
            itp = length(pfnames)-1;
            elxstr = ['/opt/elastix/bin/elastix',...
              ' -f ',p.Results.f,...
              ' -m ',p.Results.m,...
              sprintf(' -p %s',pfnames{:}),...
              ' -out ',p.Results.out];
            if ~isempty(p.Results.fMask)
                elxstr = [elxstr,' -fMask ',p.Results.fMask];
            end
            if ~isempty(p.Results.mMask)
                elxstr = [elxstr,' -mMask ',p.Results.mMask];
            end
            if ~isempty(p.Results.t0)
                elxstr = [elxstr,' -t0 ',p.Results.t0];
            end
            if ~isempty(p.Results.threads)
                elxstr = [elxstr,' -threads ',p.Results.threads];
            end
            elxstr = [elxstr,'; '];
            
            % Add transformix if Jacobian is desired 
            %   or original moving image is input:
            if ~isempty(p.Results.orig) || p.Results.jac
                if ~isempty(p.Results.orig)
                    txfn = p.Results.orig;
                else
                    txfn = p.Results.m;
                end
                elxstr = [elxstr,'/opt/elastix/bin/transformix',...
                    ' -out ',p.Results.out,...
                    ' -in ',txfn,...
                    ' -tp ',fullfile(outdir,['TransformParameters.',num2str(itp),'.txt'])];
                if p.Results.jac
                    elxstr = [elxstr,' -jac all'];
                end
                elxstr = [elxstr,'; '];
            end
            
            % Clean up temporary files
            if p.Results.CleanUp
                elxstr = [elxstr,'find ',p.Results.out,' -name "elxtemp-*" -exec rm -f {} \; '];
            end
            
            % Add command to queue file:
            qstartchk = ~exist(self.qfile,'file');
            fid = fopen(self.qfile,'a');
            if fid~=-1
                fprintf(fid,'%s\n',elxstr);
                fclose(fid);
            else
                qstartchk = false;
            end
            if qstartchk
                sysstr = ['matlab -nodesktop -nosplash -minimize '...
                        '-r "cd(''',fullfile(fileparts(which('cmi')),'elastix'),''');',...
                            'runSysQueue(''',self.qfile,''')"&'];
                system(sysstr);
            end
        end
        function deleteQFile(self)
            if exist(self.qfile,'file')
                delete(self.qfile);
            end
        end
    end
end