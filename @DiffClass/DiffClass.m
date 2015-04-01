classdef DiffClass < FitClass
    %DiffClass Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private, GetAccess=public)
        
        % ydata : Signal Intensity
        % xdata : b-value
        
    end
    
    methods
        % Constructor
        function self = DiffClass
            self.typeStr = 'Diff';
            self.defs = struct('name',{'MonoExp',...
                                       'BiExp',...
                                       'StretchedExp'},...
                               'par0',{[11000 , 0.001],...
                                       [11000 , 0.0014 , 0.0006 , 0.37],...
                                       [11000 , 0.00096 , 0.85]},...
                               'func',{@self.monexp,...
                                       @self.biexp,...
                                       @self.stretched},...
                               'labels',{{'So','ADC'},...
                                         {'So','D1','D2','f2'},...
                                         {'So','DDC','alpha'}},...
                               'bounds',{...% lsqnonlin bounds for parameter fits
                                            % [ [1xn]lowerBound ; [1xn]upperBound ]
                                         [0 0          ; 100000 0.0025],...
                                         [0 0 0 0      ; 100000 0.0025 0.0025 1],...
                                         [0 0 0        ; 100000 0.0025 1]});
            self.mind = 1;
        end
        % ~~~~ DIFFUSION MODELS --> S(b)
        function y = monexp(self,par)
            % par = [So ADC]
            S0 = par(1); ADC = par(2);
            y = S0*exp(-self.x*ADC);
        end
        function y = biexp(self,par)
            % par = [So D1 D2 f2] 
            S0 = par(1); D1 = par(2); D2 = par(3); f2 = par(4);
            y = S0*((1-f2)*exp(-self.x*D1) + f2*exp(-self.x*D2));
        end
        function y = stretched(self,par)
            % par = [So DDC alpha]
            S0 = par(1); DDC = par(2); alpha = par(3);
            y = S0*exp(-(self.x*DDC).^alpha);
        end
    end
    
end

