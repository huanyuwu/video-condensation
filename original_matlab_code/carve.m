function [X C]=carve(X,C,M,direction)
%remove one ribbon using dynamic programming
%
%input:
%X: input video (4-D array)
%C: input cost (3-D array)
%M: cumulative minimum cost
%direction: 'vert': remove vertical ribbon
%           'hor': remove horizontal ribbon
%
%output:
%X: output video (4-D array)
%C: output cost (3-D array)

H=size(X,1);
W=size(X,2);
N=size(X,4);
%C=logical(C);
%whos
switch direction
    case 'vert'
        %record the ribbon path
        %p=zeros([1,H],'single');
        [a p]=min(M(H,:,1));
        X(H,:,:,p:N-1)=X(H,:,:,p+1:N);
        C(H,:,p:N-1)=C(H,:,p+1:N);
        for i=H-1:-1:1
            p=p+M(i+1,p,2);
            X(i,:,:,p:N-1)=X(i,:,:,p+1:N);
            C(i,:,p:N-1)=C(i,:,p+1:N);
        end
        %clear M a

        %generate the video of new size (one-ribbon removal)
        %X2=zeros([H,W,3,N-1],'uint8');
        %C2=zeros([H,W,N-1],'single');
        %for i=1:H
            %X2(i,:,:,1:p(i)-1)=X(i,:,:,1:p(i)-1);
        %    X(i,:,:,p(i):N-1)=X(i,:,:,p(i)+1:N);
            %C2(i,:,1:p(i)-1)=C(i,:,1:p(i)-1);
        %    C(i,:,p(i):N-1)=C(i,:,p(i)+1:N);
        %end

    case 'hor'
        %record the ribbon path
        %p=zeros([1,W],'single');
        [a p]=min(M(W,:,1));
        X(:,W,:,p:N-1)=X(:,W,:,p+1:N);
        C(:,W,p:N-1)=C(:,W,p+1:N);
        for i=W-1:-1:1
            p=p+M(i+1,p,2);
            X(:,i,:,p:N-1)=X(:,i,:,p+1:N);
            C(:,i,p:N-1)=C(:,i,p+1:N);
        end
        %clear M a

        %generate the video of new size (one-ribbon removal)
        %X2=zeros([H,W,3,N-1],'uint8');
        %C2=zeros([H,W,N-1],'single');
        %for i=1:W
            %X2(:,i,:,1:p(i)-1)=X(:,i,:,1:p(i)-1);
        %    X(:,i,:,p(i):N-1)=X(:,i,:,p(i)+1:N);
            %C2(:,i,1:p(i)-1)=C(:,i,1:p(i)-1);
        %    C(:,i,p(i):N-1)=C(:,i,p(i)+1:N);
        %end
        
end
clear M a p direction H W N
%X=X(:,:,:,1:N-1);
%C=C(:,:,1:N-1);
X(:,:,:,end)=[];
C(:,:,end)=[];
%C=uint8(C);


