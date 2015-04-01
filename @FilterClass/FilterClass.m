classdef FilterClass < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = struct('Name',{'Gaussian','Median','Wiener','Moving Average'},...
                      'Inputs',{});
        selection
    end
    
    methods
    end
    
end

