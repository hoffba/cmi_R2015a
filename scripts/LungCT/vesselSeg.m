%vesselSeg('Z:\CT_Lung\Matlab_Scripts\Emily\BeaumontBatch\BeaumontData', 'Beaumont_sample.csv')
function [] = vesselSeg(data_dir, file_ids, resample_bool)
%vesselSeg.m Creates vessel segmentations and CSA maps for subject specified in file_ids
% Inputs:
%   data_dir : (string) Full file path to directory of folders named by subject id containing ct/lobe segmentation images
%   file_ids : (string) File path/name of csv file where first column is titled "ID" and contains all ids for which to run the pipeline
%   resample_bool : (boolean) True(default) = resample ct/lobe images
%                             False = load previously resampled files from
%                             subject directory
    
    if nargin < 2
        error('ErrorTests:paramTest', strcat("ERROR: 2 arguments must be provided to this function. \n", ...
                '1) (string) Full file path to directory of folders named by subject id containing ct/lobe segmentation images. \n', ...
                '2) (string) File path/name of csv file where first column is titled "ID" and contains all ids for which to run the pipeline. \n\n', ...
                "For example: vesselSeg('Z:/CT_Lung/Matlab_Scripts/Emily/BeaumontBatch/BeaumontData', 'Beaumont_sample.csv')"))
    elseif nargin == 2 
        resample_bool = true;
    end

    subjs = importfile(file_ids);
    for i = 1:height(subjs) 
        disp(subjs.ID(i))
        try
            subj_dir = strtrim(strcat(data_dir,  '\', string(subjs.ID(i))));
            
            %Load ct/lobes and resample if resample_bool = True
            [ct, lobes, name_sub] = getSubjectFiles(subj_dir, resample_bool);
            
            %Erode each lobe and save in subject directory
            eroded_lobes = zeros(size(lobes));
            for lobe_id = 1:5
                eroded_lobe = imerode((lobes == lobe_id), strel('sphere', 5)); 
                eroded_lobes(eroded_lobe) = lobe_id; 
            end
            eroded_lobes = uint8(eroded_lobes);
            save(strcat(subj_dir, '\', name_sub, '_erodedLobes.mat'), 'eroded_lobes');


             %Segment enhanced vessels
             vessels = vesselSeg_subj(ct, lobes);
%              vessels = vesselSeg_subj_boxes(ct, lobes); 

             %Save enhanced vessels file
            save(strcat(subj_dir, '\', name_sub, '_enhancedVessel.mat'), 'vessels');
             
             bin_vessels = binarizeVessels(vessels);
             %Save binary vessels file
             niftiwrite(uint8(bin_vessels), strcat(subj_dir, '\', name_sub, '_binVessel.nii'), 'Compressed', true);
%              save(strcat(subj_dir > 0, '\', name_sub{1}, '_binVessel_ac.mat'), 'bin_vessels');

             
             csa_map = createCSA_maps(bin_vessels);
             save(strcat(subj_dir, '\', name_sub, '_CSA_skel.mat'), 'csa_map');
             
             
             %Frangi Filter
             frangi_enhanced_vessels = vesselSeg_subj_frangi(ct, lobes); 
             save(strcat(subj_dir, '\', name_sub, '_enhancedVessel_frangi.mat'), 'frangi_enhanced_vessels');
             %bin_vessels_frangi = uint8(binarizeVessels(frangi_enhanced_vessels));
             %save(strcat(subj_dir > 0, '\', name_sub{1}, '_binVessel_frangi.nii'), 'bin_vessels_frangi');

               
             %Curvilinear Filter
             curv_enhanced_vessels = vesselSeg_subj_cruv(ct, lobes);
             save(strcat(subj_dir, '\', name_sub, '_enhancedVessel_curv.mat'), 'curv_enhanced_vessels');
             %bin_vessels_curv = uint8(binarizeVessels(curv_enhanced_vessels));
             %save(strcat(subj_dir > 0, '\', name_sub{1}, '_binVessel_curv.nii'), 'bin_vessels_curv');

        catch E
            disp(E)
        end
    end
end


function [ct, lobes, name] = getSubjectFiles(subj_dir, resample_bool)
    %Find ct/lobe filesnames
    file_ct = filterFiles(getFiles(subj_dir), 'ct', resample_bool);
    file_lobes = filterFiles(getFiles(subj_dir), 'lobe', resample_bool);

    if resample_bool
        %Resample ct/lobe images
        [ct, lobes] = resample_subj(subj_dir, file_ct, file_lobes);
        name_sub = strsplit(file_ct, '.');
        name = strcat('re_', name_sub{1});
    else
        ct = niftiread(strcat(subj_dir, '\', file_ct));
        lobes = niftiread(strcat(subj_dir, '\', file_lobes));
        name_sub = strsplit(file_ct, '.');
        name = name_sub{1};
    end
end



function[filename] = filterFiles(files, type, resample_bool)
    if strcmp(type, 'ct') && resample_bool
        files = files(contains(files, 'ct') & contains(files, 'INSP') &...
            ~contains(files, 'LD') & ~contains(files, 'COPD2') & ~contains(files, 're_'));
    
    elseif strcmp(type, 'ct')
         files = files(contains(files, 'ct') & contains(files, 'INSP') &...
            ~contains(files, 'LD') & ~contains(files, 'COPD2') & contains(files, 're_'));

    elseif strcmp(type,'lobe') 
        files = files(contains(files, 'INSP') & ~contains(files, 'LD') & ~contains(files, 'COPD2') & ...
            contains(files, 're_') & (contains(lower(files), 'lungs') | contains(lower(files), 'lobes')| contains(lower(files), 'lobe_segmentation')));
    
    elseif strcmp(type,'lobe') && resample_bool
        files = files(contains(files, 'INSP') & ~contains(files, 'LD') & ~contains(files, 'COPD2') & ...
            ~contains(files, 're_') & (contains(lower(files), 'lungs') | contains(lower(files), 'lobes')| contains(lower(files), 'lobe_segmentation')));

    else
        print('Not a valid type of file');
        assert(false);
    end
    if length(files) == 1
        filename = string(files(1));
    else
        filename = '';
    end
end
    

function [V] = vesselSeg_subj_boxes(ct, lobes)
    lobes = lobes > 0;
    I = double(ct) .*double(lobes);
    I = Normalize(I);
    clearvars ct
    
    V = zeros(size(I));
    
    %Divide I into num_boxes pieces to fit into memory
    [num_boxes, box_size] = calculateBoxSize(I);
    image_size = size(I);
    
    %MTHT3D Parameters
    no = 12; %number of orientations
    beta = 70; 
    alpha = 0.5; 
    c = 15; %parameters for the Vesselness
    alfa = -1/3; %parameters for the Neuritenees
    
    
    for i = 1:num_boxes
        if i ~= num_boxes
            box = I(:, :, 1+(i-1)*box_size(3):box_size(3)*i);
            [tmp, ~] = MTHT3D(box,0.5:1:5.5,no,beta,alpha,c,alfa);
            V(:, :, 1+(i-1)*box_size(3):box_size(3)*i) = tmp;
        else
            box = I(:, :, 1+(i-1)*box_size(3):image_size(3));
             [tmp, ~] = MTHT3D(box,0.5:1:5.5,no,beta,alpha,c,alfa);
            V(:, :, 1+(i-1)*box_size(3):image_size(3)) = tmp;
        end
    end
    V = double(V) .* double(lobes);
    
end

function [V] = vesselSeg_subj(ct, lobes)
    lobes = lobes > 0;
    I = double(ct) .*double(lobes);
    I = Normalize(I);
    clearvars ct
    
    %MTHT3D Parameters
    no = 12; %number of orientations
    beta = 70; 
    alpha = 0.5; 
    c = 15; %parameters for the Vesselness
    alfa = -1/3; %parameters for the Neuritenees
    
   [V, ~] = MTHT3D(I,0.5:1:5.5,no,beta,alpha,c,alfa); 
    V = double(V) .* double(lobes);
end

function [num_boxes, box_size] = calculateBoxSize(I)
    num_boxes = 1;
    box_size = size(I);
    while(prod(box_size) > 140000000)
        box_size(3) = floor(box_size(3)/2);
        num_boxes = num_boxes + 1;
    end
    disp(strcat('number of boxes: ', string(num_boxes)))
end

function[V_bin] = binarizeVessels(V, eroded_lobes)
%     regional_max_mask = imregionalmax(V);
%     V_bin = activecontour(V, regional_max_mask, 300, 'Chan-Vese');
    
    V_bin_orig_adp = imbinarize(V,'adaptive');
    V_bin_orig_adp_eroded = V_bin_orig_adp.*(eroded_lobes > 0);
    V_bin = activecontour(V, V_bin_orig_adp_eroded, 5, 'Chan-Vese');
end


function[V] = vesselSeg_subj_frangi(ct, lobes)
    addpath('frangi_filter_version2a')
    I = double(ct) .*double(lobes > 0);
    opt.FrangiScaleRange = [0.5 5.5];
    opt.FrangiScaleRatio = 1;
    opt.BlackWhite = false;
    [V,~,~,~,~]=FrangiFilter3D(I,opt);
    V = double(V) .* double(lobes > 0);
end


function[V] = vesselSeg_subj_cruv(ct, lobes)
    addpath('Curvilinear_Enhancement_Filter')
    I = double(ct) .*double(lobes > 0);
    V = vesselness3D(I, 0.5:1:5.5, [1,1,1], 1, true);
    V = double(V).*double(lobes > 0);
end
