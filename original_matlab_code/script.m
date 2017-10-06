eps=32;                          %stopping criterion
flexmax=3;                      %maximum value of flex parameter
colormap=zeros(2,3);            %set up colormap for output cost file
colormap(2,:)=ones(1,3);
Xinfo=aviinfo('whole_video.avi');
H=Xinfo.Height;
W=Xinfo.Width;


%flex==0
%input file: whole_video.avi, whole_cost.avi
%output file: flex0_video.avi, flex0_cost.avi
h=waitbar(0,'Computing flex=0...');
tic
X_avi=avifile('flex0_video.avi','compression','None','fps',30);
C_avi=avifile('flex0_cost.avi','colormap',colormap,'compression','None','fps',30);
for frame=1:Xinfo.NumFrames
    waitbar(frame/Xinfo.NumFrames)
    X=aviread('whole_video.avi',frame);
    X=X.cdata;
    C=aviread('whole_cost.avi',frame);
    C=C.cdata;
    %check the cost of each frame
    if sum(sum(C))>eps
    %if any(any(C))
        X_avi=addframe(X_avi,X);
        C_avi=addframe(C_avi,C);
    end
end
X_avi=close(X_avi);
C_avi=close(C_avi);
toc
close(h);

%flex>=1
%input file: flex0_video.avi, flex0_cost.avi
%output file: flex1_video.avi, flex1_cost.avi, etc.

%the following four vectors are for statistics only
vertRibbonCount=zeros(flexmax,1);
horRibbonCount=zeros(flexmax,1);
way1Count=zeros(flexmax,1);
way2Count=zeros(flexmax,1);

for flex=1:flexmax
    h=waitbar(0,['Computing flex=',num2str(flex),'...']);
    tic
    X_avi=avifile(['flex',num2str(flex),'_video.avi'],'compression','None','fps',30);
    C_avi=avifile(['flex',num2str(flex),'_cost.avi'],'colormap',colormap,'compression','None','fps',30);
    %initialization
    Mphi=flex*max(W,H)-flex+1;
    N=fix(2*Mphi);
    startframe=1;       %start reading from this frame
    endframe=N;         %stop reading until this frame
    bufferlength=0;     %the current length of the buffer
    Xinfo=aviinfo(['flex',num2str(flex-1),'_video.avi']);
    
    while startframe<=Xinfo.NumFrames
        waitbar(startframe/Xinfo.NumFrames)
        if endframe<Xinfo.NumFrames
            tempX=aviread(['flex',num2str(flex-1),'_video.avi'],startframe:endframe);
            X(:,:,:,bufferlength+1:N)=cat(4,tempX.cdata);
            clear tempX
            tempC=aviread(['flex',num2str(flex-1),'_cost.avi'],startframe:endframe);
            C(:,:,bufferlength+1:N)=cat(3,tempC.cdata);
            clear tempC
            [X C vertRibbonCount(flex) horRibbonCount(flex)]=ribboncarvemain(X,C,flex,eps,vertRibbonCount(flex),horRibbonCount(flex));        %do ribbon carving
            %[X C]=ribboncarvemainv2(X,C,flex);
            
            Np=size(X,4);
            if Np > Mphi
                Npp=Np-Mphi;
                %save the first Npp frames
                for frame=1:Npp
                    X_avi=addframe(X_avi,X(:,:,:,frame));
                    C_avi=addframe(C_avi,C(:,:,frame));
                end
                %push the remaining output to the front
                X(:,:,:,1:Mphi)=X(:,:,:,Npp+1:Np);
                C(:,:,1:Mphi)=C(:,:,Npp+1:Np);
                bufferlength=Mphi;
                %set up the next reading frames
                startframe=endframe+1;
                endframe=startframe+N-Mphi-1;
                way1Count(flex)=way1Count(flex)+1;
            else
                %no frame is saved
                bufferlength=Np;
                startframe=endframe+1;
                endframe=startframe+N-Np-1;
                way2Count(flex)=way2Count(flex)+1;
            end
            
        else
            %processing the last video chunck
            tempX=aviread(['flex',num2str(flex-1),'_video.avi'],startframe:Xinfo.NumFrames);
            X(:,:,:,bufferlength+1:bufferlength+Xinfo.NumFrames-startframe+1)=cat(4,tempX.cdata);
            clear tempX
            X=X(:,:,:,1:bufferlength+Xinfo.NumFrames-startframe+1);
            tempC=aviread(['flex',num2str(flex-1),'_cost.avi'],startframe:Xinfo.NumFrames);
            C(:,:,bufferlength+1:bufferlength+Xinfo.NumFrames-startframe+1)=cat(3,tempC.cdata);
            clear tempC
            C=C(:,:,1:bufferlength+Xinfo.NumFrames-startframe+1);
            [X C vertRibbonCount(flex) horRibbonCount(flex)]=ribboncarvemain(X,C,flex,eps,vertRibbonCount(flex),horRibbonCount(flex));
            %[X C]=ribboncarvemainv2(X,C,flex);
            
            Np=size(X,4);
            %save the entire buffer
            for frame=1:Np
                X_avi=addframe(X_avi,X(:,:,:,frame));
                C_avi=addframe(C_avi,C(:,:,frame));
            end
            clear X C
            break
        end
    end
    %the avi files of video and cost are eventually produced here
    X_avi=close(X_avi);
    C_avi=close(C_avi);
    toc
    close(h);
end
