function [img,label,fov] = readMAT(varargin)
load(varargin{1});
if ~(exist('img','var') && exist('fov','var') && exist('label','var'))
    if exist('x','var')
        if isfield(x,'idata') % Tom's cine data
            [d1,d2] = size(x(1).idata);
            [ns,nt] = size(x);
            img = zeros(d1,d2,ns,nt);
            label = cell(num2str((1:nt)'))';
            loc = zeros(ns,3);
            for i = 1:ns
                for j = 1:nt
                    img(:,:,i,j) = x(i,j).idata;
                    if (j == 1)
                        loc(i,:) = x(i,j).loc';
                    end
                end
            end
            fov = [1.25*d1 1.25*d2 ns*sqrt(sum((loc(2,:)-loc(1,:)).^2))];
        else
            img = [];
            fov = [];
            label = {};
        end
    else
        img = [];
        fov = [];
        label = {};
    end
end