% ElxClass Static function
%   Copies given text files (first column of fnames) to RegClass.odir while
%   renaming using second column of fnames. Rewrites InitialTransform 
%   parameter to link given files in new chain.
% Syntax:
%   copyTfxChain(odir,fnames)
%   copyTfxChain(odir,fnames,itname)
% Inputs:
%   odir   = (string) directory to save new Guess files into
%   fnames = (cellstr) [ copyFrom(fullpath) , newName(fullpath) ]
%   itname = (string) full filename for chain initial transform
% e.g.  RegObj0.copyTfxChain( { './TransformParameters.1.txt' , 'Guess.1.txt' ; ...
%                               './TransformParameters.0.txt' , 'Guess.0.txt' },...
%                             './InitialTransform.txt' );
function pname = copyTfxChain(odir,fnames,itname)

    if (nargin<3) || isempty(itname)
        itname = 'NoInitialTransform';
    elseif ~exist(itname,'file')
        error('File not found: %s',itname);
    end
    if ischar(fnames)
        fnames = {fnames};
    end
    nf = size(fnames,1);
    nextname = fullfile(odir,fnames{1,2});
    pname = nextname;
    for i = 1:nf
        if ~exist(fnames{i,1},'file')
            error('File not found: %s',fnames{i,1});
        else
            newname = nextname;
            if i==nf
                nextname = itname;
            else
                nextname = fullfile(odir,fnames{i+1,2});
            end
            
            % Read transform parameter file
            fid = fopen(fnames{i,1},'r');
            txt = fread(fid,'*char')';
            fclose(fid);
            
            % Change initial transform pointer
            ind = regexp(txt,'\(InitialTransformParametersFileName \"(.*?)\"\)','tokenExtents');
            txt = [txt(1:ind{1}(1)-1),nextname,txt(ind{1}(2)+1:end)];
            
            % Write new file in desired directory
            fid = fopen(newname,'w');
            fwrite(fid,txt,'char');
            fclose(fid);
            
        end
    end

end