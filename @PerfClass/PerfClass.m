classdef PerfClass < FitClass
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        defs = struct(...
              'name',{'ToftsKermode',...
                      'Patlak',...
                      'ShutterSpd',...
                      'RefRegion'},...
              'par0',{  [60 , 0.0004 , 0.004],...           % DCE - TK
                        [60 , 0.0004 , 0.004 , 0.007],...   % DCE - Patlak
                        [60 , 0.004 , 0.4 , 1],...          % DCE - ShutterSpd
                        [60 , 0.0004 , 0.004]},...          % DCE - Ref Reg
              'func',{... % function name for model
                        'tkdce','patlak','shutterspd','refreg'},...
              'labels',{{'t-delay','Ktrans','kep'},...
                        {'t-delay','Ktrans','kep','vp'},...
                        {'t-delay','L','p0','taui'},...
                        {'t-delay','Ktrans','kep'}},...
              'bounds',{...% lsqnonlin bounds for parameter fits
                           % [ [1xn]lowerBound ; [1xn]upperBound ]
                        [0 0 0      ; 90 1 1],...
                        [0 0 0 0    ; 90 1 1 1],...
                        [0 0 0 0    ; 90 1 1 5],...
                        [0 0 0      ; 90 1 1]});
                        
        % DCE Acquisition Defaults
        dceOpts = struct(...
            'r1',5.5,...        % CR relaxivity [mM^-1*s^-1]  --- BAH measured 5.5, but other refs are around 3.4
            'R1i',1/2.0,...     % T1 relaxation constant of intracellular water [s^-1]
            'R1t',1/2.2,...     % T1 relaxation constant of tissue [s^-1]
            'npools',1,...      % # of T1 pools
            'h',0.3,...         % hematocrit (in capillaries - use 0.5 for arterial blood)
            'dose',0.15,...     % mmol/kg
            'flipa',20,...      % GRE flip angle (degrees)
            'TR',0.085,...      % GRE repetition time
            'nbase',4,...       % Number of baseline images (only used to calculate baseline S)
            'naif',2,...        % AIF index (1:Maxwell; 2:Yankeelov)
            'RRktrans',0.0017,...% Ref Region Ktrans
            'RRve',0.1);        % Ref Region ve
            % Reference Region assumed kinetic values (muscle, Yankeelov 2005)
    end
    
    methods
        % Set AIF option
        function setAIF(self,val)
            if (nargin<2)
                answer = inputdlg('Input AIF index (1-2):');
                if ~isempty(answer)
                    val = str2double(answer);
                end
            end
            if isnumeric(val) && (val<3) && (val>0)
                self.naif = int8(val);
            end
        end
        % Set DCE-specific variables
        function setPerfVars(self,val)
            ok = false;
            prompt = {'r1 (/mM/s):','R1i (/s):','R1toi (/s):',...
                'Flip Angle (deg):','TR (s):','# Baseline Images:'};
            if (nargin==1) % ask user to input values
                default = {num2str(self.r1),num2str(self.R1i),num2str(self.R1t),...
                    num2str(self.flipa),num2str(self.TR),num2str(self.nbase)};
                val = str2double(inputdlg(prompt,'',1,default));
                if ~any(isnan(val))
                    ok = true;
                end
            elseif (length(val)==length(prompt)) && isnumeric(val)
                ok = true;
            end
            if ok
                self.r1 = val(1);
                self.R1i = val(2);
                self.R1t = val(3);
                self.flipa = val(4);
                self.cf = cosd(val(4));
                self.sf = sind(val(4));
                self.TR = val(5);
                self.nbase = val(6);
            end
        end
        % Set reference region data (for DCE - RefReg model)
        function setRRdata(self,y)
            if (nargin==2)
                answer = str2double(inputdlg({'Ktrans_RR:','ve_RR:'},...
                                'Reference Region Characteristics',1,...
                                {num2str(self.RRktrans),num2str(self.RRve)}));
                if ~any(isnan(answer) | (answer<0))
                    self.RRdata = y;
                    self.gdcheck(2) = false; % always load SI data, not [Gd]
                end
            end
        end
        % ~~~~ DCE-PERFUSION MODELS --> C(t)
        function y = aif(self,td)
            % td: time delay from start of imaging to bolus injection
            % returns Cp(t)
            tx = self.x-td;
            mydose = 0.15;
            switch self.naif
                case 1 % (Maxwell et al., 2002)
                    tf = mydose/0.1;
                    M = [0.8259 , 0.223 , 0.1565]*tf; % mM
                    a = [1.22 , 0.156 , 0.017]/60; % /s
                    y = M(1)*exp(-a(1)*tx) + M(2)*exp(-a(2)*tx) + M(3)*exp(-a(3)*tx);
                case 2  % Yankeelov 2003
                    tf = mydose/0.17;
                    A = 0.012*tf; B = 1.9; C = 0.05; D = 3.5*tf; E = 0.1; F = 0.0125; G = 0.84*tf;
                    y = A*(tx.^B).*exp(-C*tx) + D * (1 - exp(-E*tx)) .* exp(-F*tx) + G;
            end
            y(self.x<td) = 0;
        end
        function y = tkdce(self,par)
            % par = [td Ktrans kep]
            % OUTPUT: y = R1
            td = par(1);
            ktrans = par(2); kep = par(3);
            y = ktrans*conv(self.aif(td),exp(-kep*self.x));
            y = y(1:length(self.x));
            y = self.R1t + self.r1*y;
            y = self.sf * (1 - exp(-self.TR*y)) ./ (1 - self.cf*exp(-self.TR*y));
        end
        function y = patlak(self,par)
            % par = [td Ktrans kep vp]
            % OUTPUT: y = R1
            td = par(1);
            ktrans = par(2); kep = par(3); vp = par(4);
            taif = self.aif(td);
            y = ktrans*conv(taif,exp(-kep*self.x)); y = y(1:length(self.x));
            y = y + vp*taif;
            y = self.R1t + self.r1*y;
            y = self.sf * (1 - exp(-self.TR*y)) ./ (1 - self.cf*exp(-self.TR*y));
        end
        function y = shutterspd(self,par)
            % par = [td L p0 taui]
            %       L = kep;
            %       p0 = ve / fw --> fw = 0.8 (Yankeelov 2003)
            % OUTPUT: y = S/So
            td = par(1);
            L = par(2); p0 = par(3); taui = par(4);
            
            taif = self.aif(td);
            Co = L*conv(taif,exp(-L*self.x)); Co = Co(1:length(taif));
            g = (self.R1t - self.R1i + 1/taui)/p0;
            den = sqrt((2/taui - self.r1*Co - g).^2 + 4*(1 - p0)/(taui^2*p0));
            EL = exp(-self.TR*(2*self.R1i + self.r1*Co + g - den)/2);
            if self.npools>1
                ES = exp(-self.TR*(2*self.R1i + self.r1*Co + g + den)/2);
                aS = (1 - ( ( (self.R1i - self.R1t)/p0 - Co)*(2*p0 - 1) + 1/(taui*p0))./den)/2;
                y = self.sf  * ( (1-aS).*(1 - EL)./(1 - EL*self.cf) + aS.*(1 - ES)./(1 - ES*self.cf) );
            else
                y = self.sf * (1 - EL) ./ (1 - EL*self.cf);
            end
        end
        function y = refreg(self,par)
            % par = [td Ktrans kep]
            % OUTPUT: y = S
            td = par(1);
            ktrans = par(2); kep = par(3);
            ktransRR = self.RRktrans;
            veRR = self.RRve;
            if isempty(self.RRdata) % Estimate RR curve
                yRR = self.tkdce([td ktransRR,ktransRR/veRR]); % Using TK model (*may want to change to Patlak model in future)
            else                    % Use the given RR data curve
                yRR = self.RRdata;
            end
            R=ktrans/ktransRR;
            y = conv(yRR,exp(-kep*self.x)); y = y(1:length(self.x));
            y = R*(yRR + (ktransRR/veRR - kep)*y);
        end
    end
    
end

