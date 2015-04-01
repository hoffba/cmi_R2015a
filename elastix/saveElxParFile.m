% Function to save Elastix parameters to .txt file
function fname = saveElxParFile(s,fname)

if (nargin<2)
    [fname,path] = uiputfile('*.txt','Save Elastix Parameters');
    if fname~=0
        fname = fullfile(path,fname);
    end
end
rqflds = {'Registration','Pyramid','Metric','Optimizer','ImageSampler',...
          'InterpolatorResampler','Transform'};
if nargin && isstruct(s) && all(isfield(s,rqflds))
    fid = fopen(fname,'w');
    sect = fieldnames(s);
    for isect = 1:length(sect)
        fprintf(fid,'// ******** %s\n',sect{isect});
        vname = fieldnames(s.(sect{isect}));
        for ivname = 1:length(vname);
            FPrintPar(fid,vname{ivname},s.(sect{isect}).(vname{ivname}));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function FPrintPar(fid,name,val)
if isnumeric(val)
    vstr = ' %.3f';
    if all(val==round(val))
        vstr = ' %u';
    end
    fprintf(fid,['(%s',repmat(vstr,1,length(val)),')\n'],name,val);
elseif ischar(val)
    fprintf(fid,'(%s "%s")\n',name,val);
elseif iscellstr(val)
    fprintf(fid,['(%s',repmat(' "%s"',1,length(val)),')\n'],name,val{:});
end

