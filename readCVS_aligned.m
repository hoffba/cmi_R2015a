function PRMv3_all=readCVS_aligned






[fn,path]=uigetfile('*.csv','Select a CSV.');

fileID=fopen(fullfile(path,fn));
M = textscan(fileID,'%s %s %s','HeaderLines',1,'Delimiter',',');
fclose(fileID);

tic
count=1; % skip the header
pt_count=1;

while count<=length(M{3})
    try
        if ~strcmp(M{2}(count),'9')
            parfor j=1:3
                M{3}(count+j-1)
                data=cmi_load(1,[],M{3}(count+j-1));
                if strcmp(M{2}(count+j-1),'1')
                    data=flipdim(data,3);
                end
                if j==3
                    data(data>0.9&data<2.1)=1;data(data>1)=0;data=logical(data);
                end
                All_data(j)={data};
            end
            
            figure(100)
            for i=1:3
                subplot(1,3,i);imagesc(All_data{i}(:,:,size(All_data{i},3)-100));colormap('gray');
                ax=gca; ax.Title.String=cell2mat(M{1}(count));
            end
            pause(2)
            
            PRMv3=JDH_bounds_v3(All_data);
        else
            PRMv3=zeros(1,69);
        end
    catch err
        'Error'        
    end
    PRMv3_all(pt_count,:)=cat(2,M{1}(count),PRMv3);
    assignin('base','PRMv3_all',PRMv3_all);pause(2);
    count=count+3;
    pt_count=pt_count+1;
    
end
toc


