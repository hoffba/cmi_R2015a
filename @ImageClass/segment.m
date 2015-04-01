% ImageClass function
% Automated segmentation
function segment(self)
if self.check
    % good starting values for CT lung: seg radius=[5 5 5], Iter=100, Neighbors=0.4,
    % bottom=-100000, top=2100.
    prompt  = {'Downsize for Faster Segm','Seg Radius: single value','Stop Criteria: Iteration','Stop Criteria: % of Neighbors',...
        'Stop Criteria: Bottom Threshold','Stop Criteria: Top Threshold','Use Seed Points: 1=yes 2=no','Display Growth (y/n)'};
    title1   = 'Segmentation/Region Growing';
    lines   = (1);
    if size(handles.segm_stop,1)==2
        handles.segm_stop(3,1)=handles.thres(1);handles.segm_stop(4,1)=handles.thres(2);
    end
    % if ~isfield(handles,'segm_matrix')
        handles.segm_matrix=[size(handles.data,1) size(handles.data,2) size(handles.data,3)];
    % end
    def     = {num2str(handles.segm_matrix),...
        num2str(handles.segm_rad),num2str(handles.segm_stop(1,:)),num2str(handles.segm_stop(2,:)),...
        num2str(handles.segm_stop(3,:)),num2str(handles.segm_stop(4,:)),num2str(handles.segm_seed),handles.segm_disp};
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    input_cell = inputdlg(prompt,title1,lines,def,options);
    if isempty(input_cell)~=1
        handles.segm_stop=[];
        segm_matrix = str2num(input_cell{1}); % Slice
        if numel(segm_matrix)<3
            handles.segm_matrix(:)=segm_matrix*ones(3,1);
        else
            handles.segm_matrix=segm_matrix;
        end
        handles.segm_rad = str2num(input_cell{2}); % Slice
        handles.segm_stop(1,:) = str2num(input_cell{3}); % iterations
        handles.segm_stop(2,:) = str2num(input_cell{4}); % neighbors
        handles.segm_stop(3,:) = str2num(input_cell{5}); % lower threshold
        handles.segm_stop(4,:) = str2num(input_cell{6}); % upper threshold
        handles.segm_seed = str2num(input_cell{7}); % seed points option
        handles.segm_disp = input_cell{8};

        if numel(handles.segm_stop(1,:))>numel(handles.segm_stop(2,:))
            handles.segm_stop(2:4,:)=handles.segm_stop(2:4,end)*ones(size(handles.segm_stop(1,:)));
        elseif numel(handles.segm_stop(1,:))<numel(handles.segm_stop(2,:))
            handles.segm_stop(1,:)=handles.segm_stop(1,end)*ones(size(handles.segm_stop(2,:)));
        end
        if isempty(handles.segm_matrix)|numel(handles.segm_matrix)~=3
            handles.segm_matrix=[(size(handles.data,1)) (size(handles.data,2))...
                (size(handles.data,3))];
        end

        % interpolate Data
        if [(size(handles.data,1)) (size(handles.data,2))...
                (size(handles.data,3))]~=handles.segm_matrix
            'Interpolate Data'
            eye_orig=[(size(handles.data,1)) (size(handles.data,2))...
                (size(handles.data,3))]
            eye_ratio=(eye_orig-1)./(handles.segm_matrix-1);
            %         eye_ratio=1-eye_ratio;
            [xi,yi,zi] = meshgrid(1:eye_ratio(1,2):eye_orig(1,2)...
                ,1:eye_ratio(1,1):eye_orig(1,1)...
                ,1:eye_ratio(1,3):eye_orig(1,3));
            data= interp3(handles.data(:,:,:,handles.slider4),xi,yi,zi,'linear');
            BW1= (interp3(handles.ROI.BW(:,:,:,1),xi,yi,zi,'linear'));
            BW=BW1;BW(BW1>=0.1)=1;BW(BW1<0.1)=0;BW=logical(BW);clear BW1
            slice=round(handles.slider.slice*(handles.segm_matrix(3)/eye_orig(3)))
        else
            data=handles.data(:,:,:,handles.slider4);
            BW=logical(handles.ROI.BW(:,:,:,1));
            slice=handles.slider.slice;
            eye_orig=[(size(handles.data,1)) (size(handles.data,2))...
                (size(handles.data,3))];
        end
        if handles.segm_seed==1
            'Place seed points'
            [seed_cols seed_rows]=getpts(figure(1));
            for seed=1:length(seed_cols)
                rows=round(seed_rows(seed)*(handles.segm_matrix(1)/eye_orig(1)))
                cols=round(seed_cols(seed)*(handles.segm_matrix(2)/eye_orig(2)))
                BW(rows,cols,slice,1)=1;
            end

        end
        %----------------------------

        BWtemp=BW;
        for i=1:numel(handles.segm_stop(1,:))

            BWedge=(zeros(size(data(:,:,:))));
            BWedge=data(:,:,:);
            BWedge(data(:,:,:)>handles.segm_stop(4,i)&BWtemp==0)=0;
            BWedge(data(:,:,:)<handles.segm_stop(3,i)&BWtemp==0)=0;
            if i==1
                BWedge(BWedge>0)=1;
            else
                BWedge(BWedge>0|BWtemp==1)=1;
            end
            BWedge=logical(BWedge);
            handles.slider4_temp(1,1)=100;
            tic;
            'start growing'

            flag=0;
            roiS=0;
            count=1
            % interpolate images

            while flag==0&count<handles.segm_stop(1,i)
                BW=imdilate(BW.*handles.voi_omit_segm,ones([handles.segm_rad(i),handles.segm_rad(i),handles.segm_rad(i)]));
                BW(BWedge==0&BWtemp==0)=0;
                h = fspecial('average',[3,3]);
                BW(imfilter(BW,h)<handles.segm_stop(2,i)&BWtemp==0)=0;
                BW(imfilter(BW,h)>=0.8&BWtemp==0)=1;
                if sum(BW(BW==1))==roiS
                    flag=1;
                else
                    BWtemp=BW;
                    roiS=sum(BW(BW==1));
                    count=count+1

                end

                if strcmp(handles.segm_disp,'y')
                    BWe=edge(BW(:,:,slice,1));
                    [rows,cols]=find(BWe);
                    figure(99);imshow(data(:,:,slice)/(max(data(data>0))),'InitialMagnification','fit');colormap(handles.selected_colormap);
                    hold on
                    plot((cols), (rows),handles.voi_cs{1},'MarkerSize',handles.voi_cs{2});
                    hold off
                    set(gcf,'Name',num2str([count handles.segm_stop(1,i) handles.segm_stop(2,i)]));
                    pause(0.5)
                end
            end
        end
        if [(size(handles.data,1)) (size(handles.data,2))...
                (size(handles.data,3))]~=handles.segm_matrix
            'Interpolate Back'

            eye_ratio=(size(data)-1)./[eye_orig-1];
            %         eye_ratio=1-eye_ratio;
            [xi,yi,zi] = meshgrid(1:eye_ratio(1,2):size(data,2)...
                ,1:eye_ratio(1,1):size(data,1)...
                ,1:eye_ratio(1,3):size(data,3));
            BW1= (interp3(BW(:,:,:,1),xi,yi,zi,'nearest'));
        else
            BW1=BW;
        end
        'calculate borders'
        se90 = strel('line', 10, 90);
        se0 = strel('line', 10, 0);
        % pad to avoid image edge effects (BAH)
        [d(1) d(2) d(3)] = size(BW1);
        padl = handles.segm_rad*3;
        tmat = zeros(d+padl*2);
        tmat(padl+1:end-padl,padl+1:end-padl,padl+1:end-padl) = BW1;
        BW1 = tmat;
        BW1 = imerode(imfill(imdilate(BW1,[se90 se0]), 'holes'),[se90 se0]);
        BW1 = BW1(padl+1:end-padl,padl+1:end-padl,padl+1:end-padl);
        % undo padding
        for i=1:size(BW1,3)
            BWe=edge(BW1(:,:,i,1));
            [rows,cols]=find(BWe);
            handles.ROI.rows{i}=rows;
            handles.ROI.cols{i}=cols;
        end
        handles.ROI.BW(:,:,:,1)=BW1;size(BW1)
        guidata(gcbo,handles);
        popupmenu1_Callback;
        ['done:  ',num2str(toc/60),'  minutes']
        % ~~~~ FILL HOLES (optional, BAH) ~~~~~~~~~~~~~~~~~~~~~~~~~
        for i = 1:size(handles.ROI.BW,3)
            handles.ROI.BW(:,:,i,1) = imfill(handles.ROI.BW(:,:,i,1),'holes');
        end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
end