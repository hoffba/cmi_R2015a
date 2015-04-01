% Performs a least-squares fit on DCE-MRI timecourse data
function [pfit taif gof] = cmi_perf(t,y,model,p0,options,Crr)
% Inputs:
%           t: time (in minutes)
%           y: [CR]_toi (mM) - must be same length as t
%           model: integer flag for DCE model selection (if 0, ask for user input)
%               (1) Tofts-Kermode (TK)
%               (2) Patlak (efflux-corrected)
%               (3) Shutter-Speed (SSM)
%               (4) Reference Region - aif is the reference tissue Cref(t)
%           p0: initial parameter guess (last parameter is the taif guess)
%           Crr: (Reference Region) tissue concentration (mM) - same length as t
% Outputs:
%           pfit: vector of resulting best-fit parameters
%           taif: time delay from start of imaging to aif arrival
%           gof: goodness of fit

if nargin>=4
    % Make sure inputs are in the correct format
    t = t(:); nt = length(t);
    y = y(:);
    model = int(model);
    p0 = p0(:);
    if nargin<6
        Crr = [];
        if nargin<5
            % Use default fit options
            options = optimset( 'MaxIter',10000,...
                                'MaxFunEvals',10000,...
                                'Display','off',... % 'off' , 'on' , 'final'
                                'TolX',1e-20,...
                                'TolFun',1e-20);
        end
    end
    if length(y)==nt
        % Check for valid model selection
        if model<1
            % User-input for model
            opts = {'Tofts-Kermode','Patlak','ShutterSpeed','Ref Region'};
            [selection, ok] = listdlg('ListString',opts,'PromptString','Perfusion Model:');
        else
            ok = true;
        end
        if ok
            % Select variables based on DCE model
            switch selection
                % ~~~~~ DCE-PERFUSION MODELS ~~~~~
                case 1 % Tofts-Kermode
                    % pars = {'Ktrans','kep'};
                    % units: {1/min , 1/min}
                    func = @tkdce;
                    ub = [1 , 2];
                    lb = [0 , 0];
                    np = 2;
                case 2 % Patlak
                    % pars = {'Ktrans','kep','vp};
                    % units: {1/min , 1/min , none}
                    func = @patlak;
                    ub = [1 , 2 , 0.1];
                    lb = [0 , 0 , 0];
                    np = 3;
                case 3 % ShutterSpeed
                    % pars = {'L' , 'p0' , 'taui'};
                    % units: {1/min , none , min}
                    func = @shutterspd;
                    ub = [1 , 1, 1];
                    lb = [0 , 0 , 0];
                    np = 3;
                case 4 % Reference Region
                    % pars = {'Ktrans','kep'};
                    % units: {1/min , 1/min}
                    func = @refreg;
                    ub = [1 , 4];
                    lb = [0 , 0];
                    np = 2;
            end
            % Make sure all inputs are valid before beginning fit
            if (length(p0)==np) && ((selection~=4) || ((selection==4) && (length(Crr)==nt)))
                pfit = lsqnonlin(@(p,t0)perffitfun(p,func,t,t0,y,Crr),[p0 1],lb,ub,options);
                % Now determine the goodness of fit
                if selection==4
                    aif = Crr;
                else
                    aif = mk_aif(t,pfit(end));
                end
                yfit = feval(func,t,p,aif);
                gof = gfit(y,yfit,'4'); % 4 is RMSE
                % Now organize the outputs
                taif = pfit(end); pfit = pfit(1:end-1);
            end
        end
    end
end

% Fitting function for lsqnonlin
function F = perffitfun(p,func,t,t0,y,aif)
% Inputs:
%           p: fitting parameters
%           func: model function reference
%           t: time (min.)
%           y: timecourse data (mM)
%           aif: arterial input function (mM)
if isempty(aif)
    aif = mk_aif(t,t0);
end
F = feval(func,t,p,aif);
F = (F-y).^2;

% ~~~~ DCE-PERFUSION MODELS --> C(t)
% Arterial Input Function
function y = mk_aif(t,t0)
% returns Cp(t) -- (Maxwell et al., 2002)
% t0 refers to the time delay from imaging start to aif start
y = 0*t;
y(t>=t0) = 0.8259*exp(-1.22*(t-t0)) + 0.223*exp(-0.156*(t-t0)) + 0.1565*exp(0.017*(t-t0));
% Tofts-Kermode
function y = tkdce(t,par,aif)
% par = [Ktrans kep]
ktrans = par(1); kep = par(2);
y = ktrans*conv(aif(t),exp(-kep*t));
y = y(1:length(t));
% Efflux-Corrected Patlak
function y = patlak(t,par,aif)
% par = [Ktrans kep vp]
ktrans = par(1); kep = par(2); vp = par(3);
taif = aif(t);
y = ktrans*conv(taif,exp(-kep*t)); y = y(1:length(t));
y = y + vp*taif;
% Shutter Speed 
function y = shutterspd(t,par,aif)
% par = [L p0 taui]
L = par(1); p0 = par(2); taui = par(3);
R1toi = 60/2.2; % (min^-1) 9L tumor T1 was measured as ~2.2s (BH)
r1 = 60*5.5; % (mM^-1*min^-1) measured at 7T using calibration phantoms (BH)
R1o = 60/0.38; % hardcode for typical blood plasma relaxation in 7T
r1o = 60*3.8; % Yankeelov2003, p. 1156 (Table 1)
taif = aif(t);
g = (R1toi - R1o - 1/taui)/p0;
y = r1o*L*conv(taif,exp(-L*t)); y = y(1:length(t));
y = (2*R1toi + y - g - sqrt((2/taui + g - y).^2 + 4*(1-p0)/(taui^2*p0)) ...
    -2*R1toi)/(2*r1);
% Reference Region (Yankeelov 2005)
function y = refreg(t,par,Crr)
% par = [Ktrans kep]
ktrans = par(1); kep = par(2);
% Reference-Region parameters hardcoded
% from Yankeelov 2005 - Typical muscle values
ktransRR = 0.1; % min^-1
veRR = 0.1;
R=ktrans/ktransRR;
y = conv(Crr,exp(-kep*t)); y = y(1:length(t));
y = R*(Crr + (ktransRR/veRR - kep)*y);