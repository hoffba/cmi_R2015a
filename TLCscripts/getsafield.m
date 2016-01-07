function mf = getsafield(s,f,mtcheck)
% 'getsafield' is a Matlab10 script to extract field 'f' from structure-array 's'.
% without loops (to save time).
% Currently-designed intent is to extract numerical and character arrays from structs so
% old-style mlab scripts (based on array mathematics) can be used with
% little rewrite.  Main use is to extract image data from structures
% created by readdicom4, thus 'idata' will be default for 'f'.  eg:
% x = readdicom4; % Note, x is a structure.
% xx = getsafield(x,'idata'); Or equivalent default: xx = getsafield(x);
% returns xx as a float array with dimensions like (xres,yres,nslc,ntps).
%
% xx = getsafield(x,'loc'); % Is another example.
%
% TL Chenevert UMICH
% Copyright 2005.  The Regents of the University of Michigan.
% rev 1 TLC     August 16, 2005
% rev 2 TLC     Dec 20, 2007: Add nbyte input arg
%               nbyte = 8 = same as argument not provided, Is to return double precision floating point (matlab default)
%               nbyte = 4, Is to return single precision floating point
%               nbyte = 2, Is to return signed integer
%               You need to make calling-script properly handle requested
%               data type.
%
% TLC UMich Apr30, 2012: For some reason, a few recent PHYSIOLOGS have
% returned a few "empty cells" in the waveforms that crash the script.
% FYI, when one wave (eg. ppu) has an empty cell, they all tend to have
% that same cell index points as empty.  As a workaround, lines **80-96**
% were and "check4empty" flag was added, although the default will be to
% NOT use the workaround in order to maintain consistency with all prior use of this
% script.  In the event, if you get "CAT" errors in cell2mat, try setting
% the "check4empty" to 'y'; rerun it; then return check4empty to 'n'.
% TLC 20140409: Create an optional third input argument mtcheck = 'n' (def), or 'y' filter-out empty cells. 



subsamp = 'n'; % Normally = 'n', but n some cases, may need to subsample due to memory limits

% 20140409 Start ...
if nargin < 3
    check4empty = 'n'; % the default
else
    check4empty = mtcheck;
end

% ************************************************************************
% Unblot next line for manual override:
% check4empty = 'y'; % Default is 'n'.  Switch to 'y' if sporatic empty cells in physiolog discovered.
% 20140409 ... End
% ************************************************************************

if nargin < 2
    f = 'idata'; % Default use is to extract images after readdicom(4+)
end

struc_chk = isstruct(s); dims = size(s);

if (struc_chk == 0 || length(dims) < 2)
        mf = s;
        disp('Sorry, your input is not a structure-array. ');
        return
else
        sca = struct2cell(s);
    %     disp('FYI, this structure has fields ... ');
    %     s
end

flds = fieldnames(s);
fldi = find(strcmp(f,flds) == 1); % required field index

if (isempty(fldi))
    disp('Sorry, required field cannot be found in this structure. ');
    mf = s;
        return
    else
    %     disp('FYI, this structure has fields ... ');
    %     s
end

%if (length(dims) > 3)
 %   scaf = sca(fldi,:,:,:,:); 
if (length(dims) > 2)
    scaf = sca(fldi,:,:,:);
else
    scaf = sca(fldi,:,:);
end;
%whos scaf
%scaf{1,1,1}
    % cell array for required field
    if iscell(scaf{1,1,1}) % for the case of 1D-cell (e.g., "mark" in logvalues)
        mf = squeeze(scaf);
        mf(1:length(mf)) = mf{1:length(mf)};
    else
        
        % ********************************************
        if (check4empty == 'y')
            % For some unknown reason, a few cells may turn-up empty.  Rather
            % than crash the script, notify user and assign to "zero" or the
            % last non-empty cell value.
            if (isempty(scaf{1,1,1}))
                disp('Found empty cell in point 1');
                scaf(1,1,1) = {int16(0)};
            end % 1st point is special since no prior point available
            for iisempty = 2:length(s)
                if (isempty(scaf{1,1,iisempty}))
                    disp(['Found empty cell in point ' num2str(iisempty)]);
                    scaf(1,1,iisempty) = scaf(1,1,(iisempty-1)); 
                end % end if isempty
            end % for iisempty
        end % if check4empty
        % ********************************************
  
        
        mf = cell2mat(scaf); mf=squeeze(mf); % requested "numeric" matrix-field
        dsc = size(scaf{1,1,1});
        if ((dims(1) > 1) && (dsc(2) > 1)) % reshape for 4D-data
            mf = reshape(mf, dsc(1),dsc(2),dims(1),dims(2));
        end;
        if (length(dims) > 2)
            mf = reshape(mf, dims(3)*dims(2),1); % to comply with "getafield" 1D-output
        end;
    end;

clear s sca scaf;