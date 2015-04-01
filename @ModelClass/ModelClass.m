% Class containing methods for fitting data to a model
% DCE Models: --> All DCE data is analyzed as S/So (but input as SI)
%             (1) Tofts-Kermode 
%             (2) Efflux-corrected Patlak
%             (3) Shutter Speed
%             (4) Reference Region
% Diffusion:  (5) Mono-Exponential
%             (6) Bi-Exponential
%             (7) Stretched Exponential
classdef ModelClass < handle
    properties (SetAccess=private, GetAccess=public)
        % Model-Specific defaults & values
        defs = struct(...
            'gen',... % General curve-fitting forms
                struct(...
                    'name',{'Exponential'},...
                    'par0',{[1,1,0]},...
                    'func',{'genexp'},...
                    'labels',{{'A','B','C'}},...
                    'bounds',{[]}),...
            'perf',... % DCE-MRI models
                struct(...
                    'name',{'ToftsKermode',...
                            'Patlak',...
                            'ShutterSpd',...
                            'RefRegion'},...
                    'par0',{[30 , 0.0001 , 0.05],...           % DCE - TK
                            [30 , 0.0001 , 0.05 , 0.004],...   % DCE - Patlak
                            [30 , 0.0001 , 0.05 , 0.2],...       % DCE - ShutterSpd
                            [30 , 0.0001 , 0.05]},...           % DCE - Ref Reg
                    'func',{... % function name for model
                              'tkdce','patlak','shutterspd','refreg'},...
                    'labels',{{'t-delay','Ktrans','ve'},...
                              {'t-delay','Ktrans','ve','vp'},...
                              {'t-delay','L','p0','taui'},...
                              {'t-delay','Ktrans','ve'}},...
                    'bounds',{...% lsqnonlin bounds for parameter fits
                                 % [ [1xn]lowerBound ; [1xn]upperBound ]
                              [0 0 0      ; 200 0.1 0.5],...
                              [0 0 0 0    ; 200 0.1 0.5 0.1],...
                              [0 0 0 0    ; 200 0.1 1 5],...
                              [0 0 0      ; 200 0.1 0.5]}),...
            'diff',... % Diffusion models
                struct(...
                    'name',{'MonoExp',...
                            'BiExp',...
                            'StretchedExp',...
                            'TwoPoint'},...
                    'par0',{[11000 , 0.001],...                   % Diff - MonoExp
                            [11000 , 0.0014 , 0.0006 , 0.37],...  % Diff - BiExp
                            [11000 , 0.00096 , 0.85],...          % Diff - StretchedExp
                            [9000 , 2000 , 0.001 , 0.0002]},...   % Diff - Two Point
                    'func',{... % function name for model
                              'monexp','biexp','stretched','twopoint'},...
                    'labels',{{'So','ADC'},...
                              {'So','D1','D2','f2'},...
                              {'So','DDC','alpha'},...
                              {'Sf','Ss','ADCf','ADCs'}},...
                    'bounds',{...% lsqnonlin bounds for parameter fits
                                 % [ [1xn]lowerBound ; [1xn]upperBound ]
                              [0 0          ; 100000 0.0025],...
                              [0 0 0 0      ; 100000 0.0025 0.0025 1],...
                              [0 0 0        ; 100000 0.0025 1],...
                              nan}));
        cmodel = [1,1,1]; % index of current model for each model type
        typestrs = {'gen','perf','diff'}; % Selection for the type of model
        mtype = 1;  % index of model type
        
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
            'RRve',0.1,...       % Ref Region ve
            'R1oRR',0.7);        % Ref Region R1 (muscle)
            % Reference Region assumed kinetic values (muscle, Yankeelov 2005)
        cf  % Cosine of dceOpts.flipa
        sf  % Sine of dceOpts.flipa
        So                   % Fully relaxed magnetization signal
        normchk = false;     % Check whether ydata has been normalized to So
        
        % Data to model
        x               % [1xn] Vector of x-values (DCE: time, Diff: b-value)
        ydata           % [1xn] Vector of SI actual data
        RRdata          % [1xn] Optional vector of reference region [Gd]
        
        OptimChk = false; % Check if Optimization toolbox is available for fits
        fitOpts           % Structure of fit-function options (set in Constructor)
        % Structure of handles to plot objects
        h = struct('fig',[],'ax',[],'dline',[],'fline',[],'title',[]);
    end
    methods (Access=private)
        % Minimization function for curve fitting
        function chi = fitFunc(self,par)
            y = feval(@(x)self.(self.getDefs.func)(x),par);
            if strcmp(self.fitOpts.Display,'iter')
                self.updatePlot('fit',y);
            end
            pause(0.005)
            % For use with LSQNONLIN:
            chi = (y - self.ydata).^2;
            % For use with FMINSEARCH:
            %chi = sum(chi);
        end
        % Convert R1 values to S values
        function y = R12S(self,y,tSo)
            if nargin<3
                tSo = self.So;
            end
            y = tSo * self.sf * (1 - exp(-self.dceOpts.TR*y)) ...
                       ./ (1 - self.cf*exp(-self.dceOpts.TR*y));
        end
        % Convert S to R1 values
        function y = S2R1(self,y,tSo)
            if nargin<3
                tSo = self.So;
            end
            y = (1/self.dceOpts.TR) * log((tSo*self.sf - y*self.cf) ...
                       ./ (tSo*self.sf - y));
        end
        % Conver [Gd] to R1
        function y = gad2R1(self,y)
            y = self.dceOpts.R1t + self.dceOpts.r1*y;
        end
        % Calculate So for perfusion ydata
        function calcSo(self)
            if strcmp(self.typestrs{self.mtype},'perf') ...
                    && (length(self.ydata)>self.dceOpts.nbase)
                E = exp(-self.dceOpts.TR * self.dceOpts.R1t);
                self.So = mean(self.ydata(1:self.dceOpts.nbase)) * ...
                          (1 - E * self.cf) / ((1 - E) * self.sf);
            end
        end
        % Update plot values
        function updatePlot(self,str,y)
            go = true;
            if strcmp(str,'data')
                th = self.h.dline;
                y = self.ydata;
            elseif strcmp(str,'fit') && (nargin==3)
                th = self.h.fline;
            else
                go = false;
            end
            if go
                set(th,'YData',y);
            end
        end
    end
    methods
        % CONSTRUCTOR
        function self = ModelClass
            % Check if Optimization Toolbox is available for fits
            if license('checkout','Optimization_Toolbox')
                self.OptimChk = true;
                self.fitOpts = optimset(...
                    'MaxIter',100,...
                    'MaxFunEvals',1000,...
                    'Display','iter',... % 'off' , 'on' , 'final'
                    'TolX',1e-20,...
                    'TolFun',1e-20);
            else
                disp('Curve-Fitting not available - unable to checkout Optimization Toolbox license!')
            end
            self.sf = sind(self.dceOpts.flipa);
            self.cf = cosd(self.dceOpts.flipa);
        end
        % DESTRUCTOR
        function delete(self)
            if ishandle(self.h.fig)
                delete(self.h.fig)
            end
        end
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
        % Set current model
        function ok = setModel(self,val)
            if (nargin==1) || ~any(val==(1:length(self.cmodel)))
                % Put list of models together for user selection
                defind = self.cmodel(self.mtype);
                opts = {self.defs.(self.typestrs{self.mtype}).name};
                [val,ok] = listdlg('ListString',opts,'InitialValue',defind);
            else
                ok = true;
            end
            if ok
                self.cmodel(self.mtype) = val;
            end
        end
        % Retrieve current model type
        function str = getModType(self)
            str = self.typestrs{self.mtype};
        end
        % Retrieve current model's name
        function str = getModName(self)
            str = self.defs.(self.getModType)(self.cmodel(self.mtype)).name;
        end
        % Retrieve current model's definitions
        function str = getDefs(self)
            str = self.defs.(self.getModType)(self.cmodel(self.mtype));
        end
        % Generate current model's YData
        function y = getDefYData(self)
            tmod = self.getDefs;
            y = feval(@(x)self.(tmod.func)(x),tmod.par0);
        end
        % Retrieve current model's par0
        function val = getPar0(self)
            val = self.defs.(self.getModType)(self.cmodel(self.mtype)).par0;
        end
        % Set reference region data (for DCE - RefReg model)
        function setRRdata(self,y)
            % Input: y = S (usually muscle)
            if (nargin==2)
                % Store RR data as R1 values
                E = exp(-self.dceOpts.TR * self.dceOpts.R1oRR);
                tSo = mean(y(1:self.dceOpts.nbase)) * ...
                          (1 - E * self.cf) / ((1 - E) * self.sf);
                y = self.S2R1(y,tSo);
                self.RRdata = (y - mean(y(1:self.dceOpts.nbase))) ...
                                / self.dceOpts.r1;
            end
        end
        % Calculate goodness of fit for current model
        function gof = calcGoF(self,parin)
            if nargin==2
                y = feval(@(x)self.(self.getDefs.func)(x),parin);
                gof = sqrt(sum((y-self.ydata).^2)/length(y)); % RMSE
                gof = gof/(max(self.ydata) - min(self.ydata)); % NRMSE
            end
        end
    end
    methods % Curve calculations
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
        % ~~~~ DCE-PERFUSION MODELS --> C(t)
        function y = aif(self,td)
            % td: time delay from start of imaging to bolus injection
            % returns Cp(t)
            tx = self.x-td;
            mydose = self.dceOpts.dose;
            switch self.dceOpts.naif
                case 1 % (Maxwell et al., 2002)
                    tf = mydose/0.1;
                    M = [0.8259 , 0.223 , 0.1565]*tf; % mM
                    a = [1.22 , 0.156 , 0.017]/60; % /s
                    y = M(1)*exp(-a(1)*tx) + M(2)*exp(-a(2)*tx) + M(3)*exp(-a(3)*tx);
                case 2  % Yankeelov 2003
                    tf = mydose/0.17;
                    A = 0.012*tf;   % mM / s^B
                    B = 1.9;
                    C = 0.05;       % s^-1
                    D = 3.5*tf;     % mM
                    E = 0.1;        % s^-1
                    F = 0.0125;     % s^-1
                    G = 0.84*tf;    % mM
                    y = A*(tx.^B).*exp(-C*tx) + D * (1 - exp(-E*tx)) .* exp(-F*tx) + G;
            end
            y(self.x<td) = 0;
        end
        function y = tkdce(self,par,Cchk)
            % par = [td Ktrans ve]
            % OUTPUT: y = S
            if (nargin<3)
                Cchk = false;
            end
            td = par(1);
            ktrans = par(2); ve = par(3);
            dt = diff(self.x(1:2));
            y = ktrans*conv(self.aif(td),exp(-ktrans/ve*(self.x-td)))*dt;
            y = y(1:length(self.x));
            if ~Cchk
                y = self.R12S(self.gad2R1(y));
            end
        end
        function y = patlak(self,par)
            % par = [td Ktrans ve vp]
            % OUTPUT: y = S
            td = par(1);
            ktrans = par(2); ve = par(3); vp = par(4);
            taif = self.aif(td);
            dt = diff(self.x(1:2));
            y = ktrans*conv(taif,exp(-ktrans/ve*(self.x-td)))*dt;
            y = y(1:length(self.x)) + vp*taif;
            y = self.R12S(self.gad2R1(y));
        end
        function y = shutterspd(self,par)
            % par = [td L p0 taui]
            %       L = kep;
            %       p0 = ve / fw --> fw = 0.8 (Yankeelov 2003)
            % OUTPUT: y = S
            td = par(1);
            L = par(2); p0 = par(3); taui = par(4);
            bind = 1:self.dceOpts.nbase;
            taif = self.aif(td);
            % Calculate interstitium [Gd]
            dt = diff(self.x(1:2));
            rCo = self.dceOpts.r1 * L * conv(taif,exp(-L*(self.x-td))) * dt; rCo = rCo(1:length(taif));
            % Calculate resulting R1*
            g = (self.dceOpts.R1t - self.dceOpts.R1i + 1/taui)/p0; % s^-1
            den = sqrt( (2/taui - rCo - g).^2 + 4*(1 - p0)/(taui^2*p0) );
            y = (2*self.dceOpts.R1i + rCo + g - den)/2; % R1
            tSo = mean(self.ydata(bind));
            if self.dceOpts.npools==1
                E = exp(-self.dceOpts.TR * mean(y(bind)));
                tSo = tSo * (1-self.cf*E) / (self.sf*(1-E));
                y = self.R12S(y,tSo);
            else
                EL = exp(-self.dceOpts.TR *  y);
                ES = exp(-self.dceOpts.TR * (y + den));
                aS = (1 - ( ( (self.dceOpts.R1i - self.dceOpts.R1t)/p0 - rCo)...
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
            ktransRR = self.dceOpts.RRktrans;
            veRR = self.dceOpts.RRve;
            if isempty(self.RRdata) % Estimate RR curve
                yRR = self.tkdce([td,ktransRR,ktransRR/veRR],true);
            else
                yRR = self.RRdata;
            end
            R=ktrans/ktransRR;
            dt = diff(self.x(1:2));
            y = conv(yRR,exp(-ktrans/ve*(self.x-td)))*dt;
            y = y(1:length(self.x));
            y = R*(yRR + (ktransRR/veRR - ktrans/ve)*y);
            y = self.R12S(self.gad2R1(y));
        end
    end
end