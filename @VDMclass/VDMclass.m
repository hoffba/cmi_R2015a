classdef VDMclass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        faces    % [nf x 3] indices of vertices defining each face
        vertices % [nv x 3 x 2] locations of vertices (baseline and transformed)
        map = struct('method',{'dJ/dt','dSA/dt','GC','MC'},...
                    'label',{sprintf('$$\\dot{J}(/%s)$$',units),...
                             '$${\dot{A}(\%/yr)}$$',...
                             '$${GC}$$',...
                             '$${MC}$$'},...
                    'clim',{exp([-1,1]),[-100,100],[-0.1,0.1],[-0.3,0.3]},...
                    'logdisp',{true,false,false,false},...
                    'vals',{[],[],[],[]});
        segfname    % MHD file name of segmentation (baseline geometry)
        elxdir      % Elastix registration directory
        cfname      % TXT file name containing centerline from Mimics
        scale       % Scale factor (i.e. time, pressure, etc.)
        units       % Scale factor units
        centerLine  % Structure containing centerline measurements
    end
    events
    end
    methods
    end
end