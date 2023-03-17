function seg = lungseg_ONNX(x)
    seg = [];
    if ischar(x) && isfile(x)
        img = cmi_load(1,[],x);
    elseif isnumeric(x)
        img = x;
    else
        warning('Invalid input');
        return;
    end
    
    stat = setup_lungseg_ONNX;
    
    if stat
        d = size(img);
        img = preprocessLungCT(img,[256,256],[-1024,500]);
        seg = lungSeg(img,batchSize=8);
        seg = imresize(seg,d(1:2),'nearest');
    end
end

function stat = setup_lungseg_ONNX
    upath = userpath();
    stat = isfile(fullfile(upath,'lungmaskParams_R231.mat'));
    if ~stat
        % Need to download and set up the ONNX
        lungmask_url = "https://www.mathworks.com/supportfiles/medical/pretrainedLungmaskR231Net.onnx";
        downloadTrainedNetwork(lungmask_url,upath);
        modelfileONNX = fullfile(upath,"pretrainedLungmaskR231Net.onnx");
        modelfileM = fullfile(upath,"importedLungmaskFcn_R231.m");
        params = importONNXFunction(modelfileONNX,modelfileM);
        save(fullfile(upath,"lungmaskParams_R231.mat"),"params")

        % Need to define softmax function
        fid = fopen(modelfileM,'r+');
        str = fread(fid,'*char')';
        fclose(fid);
        str = regexprep(str,...
            '\[Vars\.x460, NumDims\.x460\] = PLACEHOLDER\(Vars\.x459\);',...
            ['Vars\.x460 = log\(softmax\(Vars\.x459,''DataFormat'',''CSSB''));\n',...
             'NumDims\.x460 = NumDims\.x459;']);
        fid = fopen(modelfileM,'w');
        fwrite(fid,str);
        fclose(fid);
        stat = true;
    end
end