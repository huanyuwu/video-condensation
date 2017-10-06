function M=cummincost(C,flex,direction)
%compute the cumulative minimum cost for all possible connected ribbons
%(using dynamic programming)
%
%input:
%C: cost array
%flex: flex parameter (non-negative integer)
%direction: 'vert': remove vertical ribbon
%           'hor': remove horizontal ribbon
%
%output:
%M(:,:,1): cumulative minimum cost (2-D array)
%M(:,:,2): direction of backtracking (2-D array)

[H W N]=size(C);
C=logical(C);

switch direction
    case 'vert'
        M=zeros([H,N,2],'single');
        sumC=permute(sum(C,2),[1 3 2]);             %sum up the cost along the 2nd dimension, and then rearrange it to 2-D array
        M(1,:,1)=sumC(1,:);
        for i=2:H
            for j=flex+1:N-flex
                [a b]=min(M(i-1,j-flex:j+flex,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-flex-1;
            end
            
            %boundary condition
            for j=1:flex
                [a b]=min(M(i-1,1:j+flex,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-j;
            end
            
            %boundary condition
            for j=N-flex+1:N
                [a b]=min(M(i-1,j-flex:N,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-flex-1;
            end
        end

    case 'hor'
        M=zeros([W,N,2],'single');
        sumC=permute(sum(C,1),[2 3 1]);         %sum up the cost along the 1st dimension, and then rearrange it to 2-D array
        M(1,:,1)=sumC(1,:);
        for i=2:W
            for j=flex+1:N-flex
                [a b]=min(M(i-1,j-flex:j+flex,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-flex-1;
            end

            %boundary condition
            for j=1:flex
                [a b]=min(M(i-1,1:j+flex,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-j;
            end

            %boundary condition
            for j=N-flex+1:N
                [a b]=min(M(i-1,j-flex:N,1));
                M(i,j,1)=sumC(i,j)+a;
                M(i,j,2)=b-flex-1;
            end
        end
        
end
