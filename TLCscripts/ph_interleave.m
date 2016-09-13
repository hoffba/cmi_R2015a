function [indx, acqord] = ph_interleave(nslice,stepsize)
% Matlab7 script to generate look-up table of the intra-TR slice interleave
% order that's implemented on Philips Achieva scanners.  This order below
% was measured emperically for multislice/multiphase 2D SShot EPI for
% perfusion imaging, TR=1.5sec, and 13-15 slices, n-nphases. For these it was
% determined that stepsize = 4.
% Thus, for a 15-slice scan, the slice order is (time going Lt to Rt):
% 1st phase=[1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12] [other phases same] ... 
% For a 14-slice scan, the slice order is:
% 1st phase=[1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 4, 8, 12] [other phases same] ...
% Purpose of this script is to properly assign sub-TR delays based on when
% each slice is acquired within each phase of the multiphase scan.  That is,  
% let delta = TR/nslice be the incremental delay from that of "Acquisition
% Time" which is the Dicom Tag used by Philips to specify time of a given
% phase.  Unfortunately, the incremental sub-TR delays dependent on
% interleave order is not available in Philips dicom tags.
% Therefore, assign td (trigger delay) value by:
% SliceLocIndx    "Trigger Delay"
%       1      AcqTime + (delta*(acqord(1) - 1))
%       2      AcqTime + (delta*(acqord(2) - 1)) 
%       3      AcqTime + (delta*(acqord(3) - 1))
%       :                :
%      indx    AcqTime + (delta*(acqord(indx) - 1))
%       :                :
%     nslice   AcqTime + (delta*(acqord(nslice) - 1))


indx = 1:nslice;
acqord = 0.*indx; % To initialize
slstart = 1;
nstep = 0;
for ii = 1:nslice
    nextslice = (stepsize*nstep) + slstart;
    if(nextslice > nslice)
        slstart = slstart+1;
        nstep = 0;
    end % end if
    nextslice = (stepsize*nstep) + slstart;
    acqord(ii) = nextslice;
    nstep = nstep+1;
    % disp(['Time Index = ' num2str(ii) ';   SliceLoc = ' num2str(acqord(ii))]);
end % end for

% for ii=1:nslice
%     disp(['Slice = ' num2str(ii) '; Delay is = ' num2str(acqord(ii)-1)]);
% end
