function [label,res] = pipeline_PTKlobes(fn)

label = [];
res = [];

try

    % Initialize 
    ptk_main = PTKMain();
    dataset = ptk_main.Load(fn);
    lobes = dataset.GetResult('PTKLobes');
    label = uint8(zeros(lobes.OriginalImageSize));
    label((1:lobes.ImageSize(1))+lobes.Origin(1)-1,...
          (1:lobes.ImageSize(2))+lobes.Origin(2)-1,...
          (1:lobes.ImageSize(3))+lobes.Origin(3)-1) = lobes.RawImage;
    label = flip(label,3);

catch err
    disp(getReport(err));
end

