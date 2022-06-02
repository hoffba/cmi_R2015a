function[V_bin] = binarizeVessels(V, eroded_lobes, ptg)
%     regional_max_mask = imregionalmax(V);
%     V_bin = activecontour(V, regional_max_mask, 300, 'Chan-Vese');
    
% Adaptive Thresholding
%     V_bin_orig_adp = imbinarize(V,'adaptive');
%
% Rosin's Thresholding
    h = imhist(V);
    T = RosinThreshold(h(2:end));
    V_bin_orig_adp = V >= (T/256);
%
% Clustering-based Thresholding
%     [~,~,LUT,~]=FastFCMeans(im2uint8(mat2gray(V)),2,4,false);
%     V_bin_orig_adp = LUT2label(im2uint8(mat2gray(V)),LUT);
%     V_bin_orig_adp_eroded = (V_bin_orig_adp == 2).*(eroded_lobes > 0);
%
% Bradley Adaptive Thresholding
%     for i = 1:size(V,3)
%         [V_bin_orig_adp(:,:,i)] = bradley(V(:,:,i),[2*floor(size(V)/16)+1], 0.5);
%     end
%    
% Bernsen Adaptive Thresholding
%     for i = 1:size(V,3)
%         [V_bin_orig_adp(:,:,i)] = bernsen(V(:,:,i),[2*floor(size(V(:,:,i))/64)+1], 0.05);
%     end
%  local contrast
% lc = 0.05;
%     V_bin_orig_adp = bernsen(V, [2*floor(size(V)/92)+1], lc);
    iter = ceil(square(ceil(ptg.*10))/2);
    if(iter > 10)
        iter = 11;
    end
    V_bin_orig_adp_eroded = V_bin_orig_adp.*(eroded_lobes > 0);
    V_bin = activecontour(V, V_bin_orig_adp_eroded, iter, 'Chan-Vese');
end