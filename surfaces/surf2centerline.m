function x = surf2centerline(f,v,c,p,n)

debug = false;

if debug
    hf = figure;
    plot3(p(:,1),p(:,2),p(:,3),'r');
    ha = gca;
    view(ha,-90,0);
    axis equal
    hold(ha,'on');
    opts = {'EdgeColor','none',...
            'FaceLighting','gouraud',...
            'FaceAlpha',0.5,...
            'AmbientStrength',0.5,...
            'DiffuseStrength',0.5,...
            'SpecularStrength',0.3,...
            'SpecularExponent',50,...
            'BackFaceLighting','reverselit'};
else
    hw = waitbar(0,'Calculating mean surface values at sectional planes ...');
end

nv = size(v,1);
np = size(p,1);
x = zeros(np,1);
valcheck = length(c)==nv;
% Loop over centerline points:
for i = 1:np
    % Find surface faces crossing normal plane:
    % a*(x-x1) + b*(y-y1) + c*(z-z1) = 0 (plane equation)
    vi = sum(bsxfun(@times,n(i,:),bsxfun(@minus,v,p(i,:))),2)>0;
    fi = find(logical(mod(sum(vi(f),2),3)));
    fcross = f(fi,:);
    uf = unique(fcross);
    
    % Find only closest connected surface:
    d = sqrt(sum(bsxfun(@minus,v,p(i,:)).^2,2));
    dd = nan(nv,1);
    dd(uf) = d(uf);
    vadd = find(dd==min(dd));
    ffi = false(size(fcross,1),1);
    while ~isempty(vadd)
        % Add faces containing vertices:
        ind = any(ismember(fcross,vadd),2);
        vadd = setdiff(fcross(ind,:),fcross(ffi,:));
        ffi(ind) = true;
    end
    fi(~ffi) = []; % Remove non-connected face indices
    
    % Display patch:
    if debug
        if i==1
            hp = patch('Faces',f(fi,:),'Vertices',v,opts{:});
        else
            hp.Faces = f(fi,:);
        end
        pause(0.1);
    else
        waitbar(i/np,hw);
    end
    
    if valcheck
        % Average vertex value:
        x(i) = mean(c(unique(f(fi,:))));
    else
        % Average face value:
        x(i) = mean(c(fi));
    end
end
if debug,delete(hf);else delete(hw),end