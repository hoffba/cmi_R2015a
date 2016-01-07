function [diffdir, diffdircosines] = getdiffgrad(info)
% This function is used to retrieve DiffusionDirectionality and direction
% cosines of diffusion gradients.  Currently only functional for Philips.
% 
% TLC 20140728.  Create script for orderly sorting of dynamic, multi-b, multi-dir unaveraged DWI.
% TLC 20141008.  Add Siemens directionality extracted from Private tags
% MD 20150706. Added GE private tags for LAB-DWI and traces
% MD 20150717. Added private tags for non-LAB DWI 
% TLC 20150827 Ancient GE data has ptag as 2-vector.
% TLC 20150827.  Switch to tryGetField(info, 'Manufacturer','UNK')
% TLC 20150904.  Add manufacturer = Imaging Biometrics LLC.

% % Initialize all manufacturer flags to 0:
% isPhilips = 0;
% isSIEMENS = 0;
% isGE = 0;
% isAgilent = 0;
% isBiometrics = 0; % 20150904
% isUnknown = 1; % 20150827
% 
% % Read manufacturer:
% % manufacturer = getdicom_field_str(info, 'Manufacturer'); %8/15/2011.
% manufacturer = tryGetField(info, 'Manufacturer','UNK'); % 20150827 
% %manufacturer = 'GE'; % 8/15/2011.
% 
% if (isempty(findstr(manufacturer ,'Philips')) == 0)
%     % Then is Philips data
%     isPhilips = 1;
%     mfg = 1;
% end
% if (isempty(findstr(manufacturer ,'SIEMENS')) == 0)
%     % Then is SEIMENS data
%     isSIEMENS = 1;
%     mfg = 2;
% end
% if (isempty(findstr(manufacturer ,'GE')) == 0)
%     % Then is GE data
%     isGE = 1;
%     mfg = 3;
% end
% if (isempty(findstr(manufacturer ,'Agilent')) == 0)
%     % Then is Agilent/Varian data
%     isAgilent = 1;
%     mfg = 10; % Skip over other future human scanners
% end
% % 201500904 start ...
% if (isempty(findstr(manufacturer ,'Biometrics')) == 0)
%     % Then is unknown data
%     isBiometrics = 1;
%     mfg = 11; % Skip over other future human scanners
% end
% % ... 20150904 end.
% % 20150827 start ...
% if (isempty(findstr(manufacturer ,'UNK')) == 0)
%     % Then is unknown data
%     isUnknown = 1;
%     mfg = 99; % Skip over other future human scanners
% end
% % ... 20150827 end.

mfg = getmfg(info);

switch mfg
    case 1 % Philips
        % Check bvalue
        bvalue = getbvalue(info);
        if (bvalue == 0)
           diffdir = 0;
           diffdircosines = [0; 0; 0];
        else
            % Set default diffdir = 1 for Isotropic, but may redefine later below:
            diffdir = 1; % Default isotropic
            diffdircosines = [0; 0; 0];
            
            % Check for "Enhanced DICOM" DiffusionDirectionality
            if (isfield(info,'DiffusionDirectionality'))
                diffdirtext = info.DiffusionDirectionality;
                if (diffdirtext(1:4) == 'NONE');
                    diffdir = 0;
                    diffdircosines = [0; 0; 0];
                elseif (diffdirtext(1:9) == 'ISOTROPIC')
                    diffdir = 1;
                    diffdircosines = [0; 0; 0];
                elseif (diffdirtext(1:11) == 'DIRECTIONAL')
                    diffdir = 2;
                    if (isfield(info,'DiffusionGradientDirectionSequence'))
                        if (isfield(info.DiffusionGradientDirectionSequence,'Item_1'))
                            if (isfield(info.DiffusionGradientDirectionSequence.Item_1,'DiffusionGradientOrientation'))
                                diffdircosines = info.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation;
                            end % if DiffusionGradientOrientation
                        end % if Item_1
                    end % if DiffusionGradientDirectionSequence
                end % if diffdirtext
                
            % OK, try "Non Enhanced DICOM" assumption  
            elseif (isfield(info,'DiffusionGradientOrientation'))
                diffdircosines = info.DiffusionGradientOrientation;
                if ( isequal(diffdircosines,[0;0;0]) )
                    diffdir = 1;
                else % then at least one diffdircosines is not zero
                    diffdir = 2;
                end % isequal
            else % Unexplained non-zero bvalue yet 'DiffusionDirectionality' or 'DiffusionGradientOrientation' do not exist
               % disp('Oops - something wrong in getdiffgrad DICOM ... ');
            end % if DiffusionDirectionality
        end % if bvalue
        
    case 2 % Siemens
        
        % Start 20141008 ...  First blot defaults diffdir and
        % diffdircosines:
%         % WIP: for now ASSUME no DiffusionDirectionality
%         diffdir = 0;
%         diffdircosines = [0; 0; 0];

        
        % Check bvalue
        bvalue = getbvalue(info);
        if (bvalue == 0)
           diffdir = 0;
           diffdircosines = [0; 0; 0];
        else
           % Set default diffdir = 1 for Isotropic, but may redefine later below:
           diffdir = 1; % Default isotropic
           diffdircosines = [0; 0; 0];
           
           % Check for "Private_0019_100e", suggesting directionality.
           % FYI: Siemens does NOT also include the trace DWI with directional ones.
           % The trace DWI is probably in the subsequent series.
           if (isfield(info,'Private_0019_100e'))
               diffdir = 2; % Then is directional
               ptag = info.Private_0019_100e;
               if(isa(ptag,'uint8'))
                   lnptag = length(ptag);
                   if (lnptag >= 24)
                       diffvector = [typecast(ptag(1:8),'double'); typecast(ptag(9:16),'double'); typecast(ptag(17:24),'double')];
                   else
                       disp('Opps - something wrong in DICOM. ... ');
                   end % if lnptag
               else % Then ptag is likely already class "double"
                   diffvector = ptag;
               end % if isa ptag
               
               diffmagn = sqrt( sum(diffvector.^2) );
               diffdircosines = diffvector/diffmagn;
               
           end % if isfield
        end % if bvalue
        % ... End 20141008.
        
    case 3 % GE - see MD 20150706  
        % WIP: for now ASSUME no DiffusionDirectionality
        diffdir = 0;
        diffdircosines = [0; 0; 0];
                
        % Check bvalue  % start MD 20150706 (assume simple DICOM)
        bvalue = getbvalue(info);
        if (bvalue == 0)
           diffdir = 0;
           diffdircosines = [0; 0; 0];
        else % b>0
           % Set default diffdir = 1 for Isotropic, but may redefine later below:
           diffdir = 1; % Default isotropic (ptag = 15)
           diffdircosines = [0; 0; 0];
           
           if (isfield(info,'Private_0043_1030'))
              % diffdir = 2; % Then is directional
               ptag = info.Private_0043_1030;
               % switch ptag  % DWI-DIR
                switch ptag(1) % 20150827 try 1st element only.  DWI-DIR
                   case 3   % 3,4,5: LAB-DWI
                       diffdir = 2;
                       diffdircosines = [1; 0; 0]; % R/L-DWI
                   case 4
                       diffdir = 2;
                       diffdircosines = [0; 1; 0]; % A/P-DWI
                   case 5
                       diffdir = 2;
                       diffdircosines = [0; 0; 1]; % S/I-DWI
                   case 16   % non-LAB DWI % MD 20150717
                       diffdir = 2;
                       if (isfield(info,'Private_0019_10bb') && isfield(info,'Private_0019_10bc') && isfield(info,'Private_0019_10bd')) 
                           diffdircosines = [info.Private_0019_10bb; info.Private_0019_10bc; info.Private_0019_10bd];
                       end
                   otherwise % isotropic
                       diffdircosines = [0; 0; 0]; % trace-DWI (case 14/15)
               end % LAB-DWI switch
%            elseif (isfield(info,'Private_0019_10bb') && isfield(info,'Private_0019_10bc') && isfield(info,'Private_0019_10bd')) 
%                % non-LAB DWI
%                diffdir = 2;
%                diffdircosines = [info.Private_0019_10bb; info.Private_0019_10bc; info.Private_0019_10bd];
            end % if isfield
        end % if bvalue
        % ... End 20150706.
        
    otherwise % TO BE ADDED LATER
        % WIP: for now ASSUME no DiffusionDirectionality
        diffdir = 0;
        diffdircosines = [0; 0; 0];
       
end % end switch mfg

return;
