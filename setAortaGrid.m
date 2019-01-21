% Script to set aorta registration grid spacing to consistent spatial lengths
function setAortaGrid(regObj)

if length(regObj.elxObj.Schedule)==2 && length(regObj.elxObj.Schedule{2}.GridSpacingSchedule)==9

    grdsp = [12,12,12,6,6,6,3,3,3]; % grid spacing, in mm
    regObj.UDpreproc('dilaten',round(repmat(grdsp(4:6),2,1) ./ [regObj.cmiObj(1).img.voxsz; regObj.cmiObj(2).img.voxsz]));
    regObj.setElxPar(2,'GridSpacingSchedule',round(grdsp./repmat(regObj.cmiObj(1).img.voxsz,1,3)));

else
    error('Invalid Elastix parameters. Requires (1) Affine and (2) Warp with 3 resolutions.')
end