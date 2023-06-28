function pipeline_save_fig(hf,fname)

if ischar(fname)
    fname = {fname};
end
for i = 1:numel(fname)
    [fpath,~,ext] = fileparts(fname{i});
    if ~isfolder(fpath)
        mkdir(fpath);
    end
    switch ext
        case {'.tif','.tiff'}
            print(hf,fname{i},'-dtiff');
        case '.gif' % Collate into GIF series
            collate_fig(hf,fname{i});
    end
end