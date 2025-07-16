function seg = pipeline_PTKlobes(fn,fn_log)

try

    % Initialize 
    ptk_main = PTKMain(CoreReporting([],false,fn_log));
    dataset = ptk_main.Load(fn);
    lobes = dataset.GetResult('PTKLobes');
    seg = uint8(zeros(lobes.OriginalImageSize));
    seg((1:lobes.ImageSize(1))+lobes.Origin(1)-1,...
          (1:lobes.ImageSize(2))+lobes.Origin(2)-1,...
          (1:lobes.ImageSize(3))+lobes.Origin(3)-1) = lobes.RawImage;
    seg = flip(seg,3);

    % Cleanup
    ptk_main.DeleteCacheForAllDatasets;

catch err
    seg = getReport(err);
end


