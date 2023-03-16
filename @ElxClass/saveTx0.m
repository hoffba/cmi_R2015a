% ElxClass function
% Save Elastix Initial Transform file
function fname = saveTx0(self,x,~)

fname = '';
if ~isempty(self.Tx0)
    if (nargin==3) && ishandle(x)
        [fname,path] = uiputfile('*.txt','Save Initial Transform',...
                                         'InitialTransform.txt');
        if fname~=0
            fname = fullfile(path,fname);
        else
            fname = '';
        end
    elseif ~ischar(x)
        fname = '';
    else
        fname = x;
    end
    if ~isempty(fname)
        permchk = strcmp(self.outfmt,'.mhd');
        switch self.Tx0.Transform
            case 'TranslationTransform'
                pari = [2  1  3];
            case 'EulerTransform'
                pari = [2  1  3  5  4  6];
            case 'SimilarityTransform'
                pari = [2  1  3  5  4  6  7];
            case 'AffineTransform'
                pari = [5  4  6  2  1  3  8  7  9  11  10  12];
        end
        fid = fopen(fname,'w');
        if fid>2
            fstr = fieldnames(self.Tx0);
            for i = 1:length(fstr)
                val = self.Tx0.(fstr{i});
                if permchk && strcmp(fstr{i},'TransformParameters')
                    val = val(pari);
                end
                vstr = '';
                if ischar(val)
                    vstr = [' "',val,'"'];
                elseif iscellstr(val)
                    vstr = sprintf(' "%s"',val{:});
                elseif all(round(val)==val)
                    vstr = sprintf(' % .0f',val);
                elseif isnumeric(val)
                    vstr = sprintf(' % .8f',val);
                end
                fprintf(fid,'(%s%s)\n',fstr{i},vstr);
            end
            fclose(fid);
        else
            error(['Could not open file: ',fname])
        end
    end
end