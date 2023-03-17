function volInp = preprocessLungCT(V,imSize,th)
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %
% Prepare Lung CT image for semantic segmentation
%
% Input:
%   V: 2-D slice or 3-D volume of a LungCT scan
%
%  imSize: two-element vector containing the target size of the first two
%       dimensions to be used for resizing the image
%
%   th: two-element array containing the threshold to be used for
%       preprocessing the image
%
% Output:
%   volInp: 2-D slice or 3-D volume (same size as input V) with x and y
%           dimensions 256x256 of single class (normalized values)
%   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %   %

% Cast Volume to Single Type
V = single(V);

% Resize the first two dimensions of V to size imSize
V = imresize(V,[imSize(1) imSize(2)]);

% Threshold and normalize intensity values
volInp = rescale(V,0,1,InputMin=th(1),InputMax=th(2));
end