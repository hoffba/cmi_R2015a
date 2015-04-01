classdef DCEclass < FitClass
    %DCEclass Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private, GetAccess=public)
        
        % ydata : Signal Intensity
        % xdata : Seconds
        
        % Acquisition parameters:
        TR = 0.085; % (sec)
        flipa = 20; % (deg)
        cf          % calculate cosd(flipa) and sind(flipa)
        sf
        
        % Assumed parameters:
        r = 5.5;       % [Gd] relaxivity @9.4T (measured by BAH)
        Ri = 1/2;      % Intracellular water relax rate (/sec), mouse tumor @9.4T
        Rt = 1/2.2;    % Tissue relax rate (/sec) for mouse tumor @9.4T
        Rb = 1/1.664;  % Blood relaxation rate (/sec) for Human @7.0T
        npools = 1;
        hct = 0.42;
        nbase = 4;
        
        % Arterial Input Function (AIF)
        aifType = 1;    % Type of AIF being used:
                        %   0 = empirical
                        %   1 = mouse1
                        %   2 = mouse2
                        %   3 = rat
                        %   4 = human
        dose = 0.115;   % (mg/kg)
        aif             % vector of Cp values ([Gd],mM)
        
        % Reference Region parameters:
        % Reference Region assumed kinetic values (muscle, Yankeelov 2005)
        RtRR = 0.7;    % Muscle relax rate (/sec) for rat @9.4T
        RRktrans = 0.0017; % (/sec)
        RRve = 0.1;
        
    end
    
    methods
        % Contructor
        function self = DCEclass
            self.cf = cosd(self.flipa);
            self.sf = sind(self.flipa);
            self.typeStr = 'DCE';
            self.defs = struct('name',{'ToftsKermode',...
                                       'Patlak',...
                                       'ShutterSpd',...
                                       'RefRegion'},...
                               'par0',{[30 , 0.0001 , 0.05],...
                                       [30 , 0.0001 , 0.05 , 0.004],...
                                       [30 , 0.0001 , 0.05 , 0.2],...
                                       [30 , 0.0001 , 0.05]},...
                               'func',{... % function name for model
                                         @self.tkdce,...
                                         @self.patlak,...
                                         @self.shutterspd,...
                                         @self.refreg},...
                               'labels',{{'t-delay','Ktrans','ve'},...
                                         {'t-delay','Ktrans','ve','vp'},...
                                         {'t-delay','L','p0','taui'},...
                                         {'t-delay','Ktrans','ve'}},...
                               'bounds',{...% lsqnonlin bounds for parameter fits
                                            % [ [1xn]lowerBound ; [1xn]upperBound ]
                                         [-100 0 0      ; 200 0.1 0.5],...
                                         [-100 0 0 0    ; 200 0.1 0.5 0.1],...
                                         [-100 0 0 0    ; 200 0.1 1 5],...
                                         [-100 0 0      ; 200 0.1 0.5]});
            self.mind = 1;
        end
        % Convert R1 values to S values
        function y = R2S(self,y,So)
            if nargin>1
                y = So * self.sf * (1 - exp(-self.TR*y)) ...
                        ./ (1 - self.cf*exp(-self.TR*y));
            end
        end
        % Convert S to R1 values
        function y = S2R(self,y,So)
            if nargin==3
                y = (1/self.TR) * log((So*self.sf - y*self.cf) ...
                           ./ (So*self.sf - y));
            end
        end
        % Convert [Gd] to R1
        function y = gad2R(self,y)
            y = self.Rt + self.r*y;
        end
        % Convert R1 to [Gd]
        function y = R2gad(self,y)
            y = (y - self.Rt)/self.r;
        end
        % Calculate So for perfusion ydata
        function So = calcSo(self,y,R)
            if (nargin<3)
                R = self.Rt;
            end
            if (nargin<2)
                y = self.ydata;
            end
            if (length(y)>self.dceOpts.nbase)
                E = exp(-self.TR * R);
                So = mean(y(1:self.nbase)) * ...
                          (1 - E * self.cf) / ((1 - E) * self.sf);
            end
        end
        % Shift AIF in time
        function y = shiftAIF(self,td)
            if self.mind
                y = self.calcAIF(td);
            else
                y = interp1(self.xdata,self.aif,self.xdata+td,'spline');
                y(isnan(y)) = 0; % Extrapolates to 0
            end
        end
        
        
        % ~~~~ DCE-PERFUSION MODELS
        function y = calcAIF(self,td)
            % Calculate AIF from model:
            tx = self.xdata - td;
            switch self.mind
                case 1 % (Maxwell et al., 2002)
                    M = [8.259 , 2.23 , 1.565]*self.dose; % mM
                    a = [1.22 , 0.156 , 0.017]/60; % /s
                    y = M(1)*exp(-a(1)*tx) + M(2)*exp(-a(2)*tx) ...
                                    + M(3)*exp(-a(3)*tx);
                case 2  % Yankeelov 2003
                    tf = self.dose/0.17;
                    A = 0.012*tf;   % mM / s^B
                    B = 1.9;
                    C = 0.05;       % s^-1
                    D = 3.5*tf;     % mM
                    E = 0.1;        % s^-1
                    F = 0.0125;     % s^-1
                    G = 0.84*tf;    % mM
                    y = A*(tx.^B).*exp(-C*tx) ...
                             + D * (1 - exp(-E*tx)) .* exp(-F*tx) + G;
                case 3
                case 4
            end
            y(self.xdata<td) = 0;
        end
        function y = tkdce(self,par)
            % par = [td Ktrans ve]
            % OUTPUT: y = S
            td = par(1);
            ktrans = par(2); ve = par(3);
            dt = diff(self.xdata(1:2));
            y = ktrans*conv(self.shiftAIF(td),...
                            exp(-ktrans/ve*(self.xdata-td)))*dt;
            y = y(1:length(self.xdata));
        end
        function y = patlak(self,par)
            % par = [td Ktrans ve vp]
            % OUTPUT: y = S
            td = par(1);
            taif = self.shiftAIF(td);
            ktrans = par(2); ve = par(3); vp = par(4);
            dt = diff(self.xdata(1:2));
            y = ktrans*conv(taif,exp(-ktrans/ve*(self.xdata-td)))*dt;
            y = y(1:length(self.xdata)) + vp*taif;
            y = self.R2S(self.gad2R(y));
        end
        function y = shutterspd(self,par)
            % par = [td L p0 taui]
            %       L = kep;
            %       p0 = ve / fw --> fw = 0.8 (Yankeelov 2003)
            % OUTPUT: y = S
            td = par(1);
            taif = self.shiftAIF(td);
            L = par(2); p0 = par(3); taui = par(4);
            bind = 1:self.dceOpts.nbase;
            % Calculate interstitium [Gd]
            dt = diff(self.x(1:2));
            rCo = self.r * L * conv(taif,exp(-L*(self.xdata-td))) * dt;
            rCo = rCo(1:length(taif));
            % Calculate resulting R1*
            g = (self.Rt - self.Ri + 1/taui)/p0; % s^-1
            den = sqrt( (2/taui - rCo - g).^2 + 4*(1 - p0)/(taui^2*p0) );
            y = (2*self.Ri + rCo + g - den)/2; % R1
            tSo = mean(self.ydata(bind));
            if self.npools==1
                E = exp(-self.TR * mean(y(bind)));
                tSo = tSo * (1-self.cf*E) / (self.sf*(1-E));
                y = self.R2S(y,tSo);
            else
                EL = exp(-self.TR *  y);
                ES = exp(-self.TR * (y + den));
                aS = (1 - ( ( (self.Ri - self.Rt)/p0 - rCo)...
                            *(2*p0 - 1) + 1/(taui*p0))./den) / 2;
                bEL = mean(EL(bind)); bES = mean(ES(bind)); baS = mean(aS(bind));
                tSo = tSo ./ ( ((1-baS)*self.sf*(1-bEL)) / (1-bEL*self.cf) ...
                             + (   baS *self.sf*(1-bES)) / (1-bES*self.cf) );
                y = tSo * self.sf  * ( (1-aS).*(1 - EL)./(1 - EL*self.cf) ...
                                             + aS.*(1 - ES)./(1 - ES*self.cf) );
            end
        end
        function y = refreg(self,par)
            % par = [td Ktrans ve]
            % OUTPUT: y = S
            td = par(1);
            ktrans = par(2); ve = par(3);
            if isempty(self.aif) % Estimate RR curve
                yRR = self.tkdce([td,ktransRR,ktransRR/veRR],true);
            else
                yRR = self.shiftAIF(td);
            end
            R=ktrans/self.RRktrans;
            dt = diff(self.x(1:2));
            y = conv(yRR,exp(-ktrans/ve*(self.x-td)))*dt;
            y = y(1:length(self.x));
            y = R*(yRR + (self.RRktrans/self.RRve - ktrans/ve)*y);
            y = self.R2S(self.gad2R(y));
        end
    end
    
end

