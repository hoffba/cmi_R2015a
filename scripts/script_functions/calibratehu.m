function calibrated_datavol = calibratehu( datavol, mean_air, mean_blood,...
                                        target_air, target_blood)
%calibratehu_ Rescale CT HU data volume 
%   
%Description:
% Take a data set with values from 0 on up and convert to actual
% HU values by subtracting -1024 and then doing a conversion
% based on the equation Vhu = Vinit*m + C where m is a slope and C
% is a constant determined from the following equations:
%   Cblood = desired value of blood in aorta (30-45 HU)
%   Cair = desired value of air above abdomen (-1024 HU GE)
%   Cmult = (Cblood - Cair)/(Mblood - Mair)
%   Cadd = Cblood - m*Mblood (or C = Cair = m*Mair)
%
% Revision history:
% 26 January 2012  jlb  Created to do a correction of CT values
% 22 March 2012    jlb  Fixed m calc because Blood & Air reversed
% 23 March 2012    jlb  Add in Targets as Parameters
% 26 March 2012    jlb  Change to ASSUME HU input
% --------------------------------------------------------------------

%
% Set constants for desired target means of blood and air:
%
% Typical expected CT HU values:
% Cblood = 40;
% Cair = -1000;
%
% COPD Gene 120v/0.625slice thickness/Body/Bone typical mean values
% Cblood=37 and Cair=-995
% 
% These are the target values for correction
%Cblood = 37;
%Cair = -995;

Cblood = target_blood;
Cair = target_air;

%
% Process input

Mblood = mean_blood;
Mair = mean_air;

% Compute Cmult and Cadd
Cmult = (Cair - Cblood)/(Mair - Mblood);
Cadd = Cblood - Cmult*Mblood;

% Cycle through all image vectors and process
for i=1:size(datavol,4)
    % Assume data in HU so DON'T Subtract GE -1024 offset first
    %datavol(:,:,:,i) = -1024 + datavol(:,:,:,i);
    
    % Do multiply by scaling and then offset addition
    datavol(:,:,:,i) = Cmult*datavol(:,:,:,i);
    datavol(:,:,:,i) = datavol(:,:,:,i) + Cadd;
    
    % Assume data in HU so DON'T Add back GE -1024 offset to make values positive
    %datavol(:,:,:,i) = 1024 + datavol(:,:,:,i);
    
end

% show that we did work!
disp('Calibrate HU using mult & add:')
disp(Cmult);
disp(Cadd);

% return corrected volume
calibrated_datavol = datavol;

end

