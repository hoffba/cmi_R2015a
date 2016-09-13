function handles=DT_IM(handles);

sq=size(handles.data.A);
h = waitbar(0,'Calculate IM');

stack=handles.stack;
if isfield(handles,'num_slices')==1
    SN=handles.num_slices;
else
    SN=handles.SN;
end
handles.data.IM=zeros(size(handles.data.A,1),size(handles.data.A,2),SN,handles.bV(1,1),handles.DWI); % preallocate
for st=1:stack
    if stack>1
        if st==1
            SNi=1;
            SNf=SN/stack
        else
            SNi=SN/stack+1;
            SNf=SN
        end
    else
        SNi=1;
        SNf=SN;
    end
    for q=1:handles.bV(1,1)
        for m=1:handles.DWI
            for n=SNi:SNf
                waitbar((st*q*m*n)/(stack*sq(3)*sq(5)*sq(4)))
                M=handles.data.A(:,:,n,q,m)./(handles.data.A1(:,:,n)+realmin);
                N=log(M+realmin);
                clear M
                handles.data.IM(:,:,n,q,m)=N;
                clear N
            end
        end
    end
end
close(h)

