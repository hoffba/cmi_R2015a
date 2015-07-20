function C = PRMcalcMF(obj)

C = {};
prm = [];
switch class(obj)
    case 'CMIclass'
        prm = obj.img.prm.mat;
        obj = obj.img.prm;
    case 'ImageClass'
        prm = obj.prm.mat;
        obj = obj.prm;
    case 'PRMclass'
        prm = obj.mat;
end

if ~isempty(prm) && obj.check
    thr = 1:max(prm(:));
    [MF,labels] = minkowskiFun(prm,thr,'==','voxsz',obj.voxsz,'prog',true);
    [prml,prmv] = obj.getStats;
    C = [ {[]} ,    {'PRM%'} ,          labels ;...
          prml(:) , num2cell(prmv(:)) , num2cell(MF)];
    
    disp(mat2clip(C));
    disp(' ... Copied to clipboard.');
end

