% ImageClass function
function [m,b] = CTphantCal(self,vec,opt)
% Must have phantom rods circled on both ends (2 end slices)

m = [];
b = [];
phants = {'QCT Pro',...
          'MicroMouse-small',...
          'MicroMouse-large',...
          'Simple Air/Water/Bone'};
ok = true;
if (nargin~=3) || isempty(opt)
    [opt,ok] = listdlg('ListString',phants);
elseif ~any(opt==(1:length(phants)))
    ok = false;
end
slc = find(any(any(self.mask.mat,1),2));
if ok 
    if (nargin==1) || ~any(vec==(1:self.dims(4)))
        vec = 1;
    end
    % Determin rod properties:
    switch opt
        case 1 % Joe's QCT Pro
            rho = [-51.83 , -53.4 , 58.88 , 157.05 , 375.83]; % mg/cc K2HPO4
            rhoW = [1012.25 , 1056.95 , 1103.57 , 1119.52 , 923.2]; % mg/cc H2O
            rnames = {'rodA','rodB','rodC','rodD','rodE'};
            r = 1.5 / self.voxsz(1); % (mm) exclude outer volume of the rod
            zoff = [-0.2174 ; 999.6];
        case 2 % MicroMouse090 small rods
            rho = [0 , 50 , 100, 250 , 500, 750]; % mg/cc Hydroxyapatite
            rhoW = 0;
            rnames = {'rod000','rod050','rod100','rod250','rod500','rod750'};
            r = 0.5 / self.voxsz(1); % (mm) exclude outer volume of the rod
            zoff = 0;
        case 3 % MicroMouse090 large rods
            rho = [0 , 50 , 250 , 750]; % mg/cc Hydroxyapatite
            rhoW = 0;
            rnames = {'rod000','rod050','rod250','rod750'};
            r = 1 / self.voxsz(1); % (mm) exclude outer volume of the rod
            zoff = 0;
        case 4 % Simple Water/Bone
            rho = [ 0 , 1000 ]; % not sure if the gammex 450 mineral density is correct
            rhoW = 0;
            rnames = {'water','G450'};
            r = 0.75 / self.voxsz(1); % (mm) exclude outer volume of the rod
            zoff = 0;
    end
    % Find centers of ROIs
    z = [slc(1),slc(end)];
    nrods = length(rho);
    tmask1 = self.mask.mat(:,:,z(1));
    tmask2 = self.mask.mat(:,:,z(2));
    cc1 = bwconncomp(tmask1);
    cc2 = bwconncomp(tmask2);
    if (cc1.NumObjects==nrods) && (cc2.NumObjects==nrods)
        % Find center and mean of each ROI
        timg1 = self.mat(:,:,z(1),vec);
        timg2 = self.mat(:,:,z(2),vec);
        x1 = zeros(1,nrods); y1 = x1; x2 = x1; y2 = x1;
        rmeans = zeros(2,nrods);
        for i = 1:nrods
            [y,x] = ind2sub(self.dims(1:2),cc1.PixelIdxList{i});
            x1(i) = mean(x); y1(i) = mean(y);
            rmeans(1,i) = mean(timg1(cc1.PixelIdxList{i}));
            [y,x] = ind2sub(self.dims(1:2),cc2.PixelIdxList{i});
            x2(i) = mean(x); y2(i) = mean(y);
            rmeans(2,i) = mean(timg2(cc2.PixelIdxList{i}));
        end
        % Sort by image value so ROIs match up
        [~,torder] = sort(rmeans(1,:)); x1 = x1(torder); y1 = y1(torder);
        [~,torder] = sort(rmeans(2,:)); x2 = x2(torder); y2 = y2(torder);
        % Generate and save cylinder VOIs to measure rod values
        [xi,yi,zi] = meshgrid(1:self.dims(2),1:self.dims(1),1:self.dims(3));
        zi = (zi - z(1)) / diff(z);
        tmask = false(self.dims(1:3));
        omask = tmask;
        rmeans = zeros(1,nrods);
        for i = 1:nrods
            x = zi * (x2(i) - x1(i)) + x1(i);
            y = zi * (y2(i) - y1(i)) + y1(i);
            tmask = ((sqrt((xi-x).^2 + (yi-y).^2)<r) & (zi>=0) & (zi<=1));
            saveMASK(fullfile(self.dir,[self.name,'_',rnames{i},'.mask']),tmask);
            omask = omask | tmask;
            tmask = find(tmask);
            rmeans(i) = mean(self.mat(tmask + (vec-1)*prod(self.dims(1:3))));
            fprintf('Mean - %6s : %8.2f\n',rnames{i},rmeans(i));
        end
        % Scale image to convert to density values
        b = (rmeans - rhoW)';
        A = [rho',ones(nrods,1)];
        z = (A'*A) \ A' * b + zoff; % [sigma,beta]_ref
        m = 1/z(1);
        b = -z(2)/z(1);
        fprintf('Image Scale (y=m*x+b): m = %6.4f , b = %6.2f\n',m,b);
        figure,plot(rmeans,rho,'ob',rmeans,rmeans*m+b,'r')
        M = m * self.scaleM(vec);
        B = m * self.scaleB(vec) + b;
        self.imgScale(vec,M,B)
        self.mask.merge('replace',omask); % set to current mask
    else
        error('# VOIs does not match # of rods!')
    end
end
