function [X C vertRibbonCount horRibbonCount]=ribboncarvemain(X,C,flex,eps,vertRibbonCount,horRibbonCount)
%The main function for ribbon carving
%
%input:
%X: input video (4-D array)
%C: input cost (3-D array)
%flex: flex parameter
%eps: stopping criterion
%
%output:
%X: output video (4-D array)
%C: output cost (3-D array)
%
%the number of carved-out vertical ribbons in this function will be
%vertRibbonCount(output) - vertRibbonCount(input).
%the number of carved-out horizontal ribbons in this function will be
%horRibbonCount(output) - horRibbonCount(input).
%
%PATH: vector representing the order of horizontal or vertical ribbon removal.
%      The elements are either 1 or 2:
%      1: remove horizontal ribbon
%      2: remove vertical ribbon

%Huan-Yu Wu


[H W rgb N]=size(X);
%PATH=zeros(1,N);
%i=1;
%X2=X;
%C=costBS_MRF(X2,theta,gamma,alpha,iteration);                %compute the cost function (only once)

My=cummincost(C,flex,'vert');      %compute the cumulative minimum cost
[Cv b]=min(My(H,:,1));                 %the cost of removing a vertical ribbon
Mx=cummincost(C,flex,'hor');
[Ch a]=min(Mx(W,:,1));                 %the cost of removing a horizontal ribbon
%Ch
%a
%Cv
%b
while Ch<=eps || Cv<=eps
    %remove a vertical ribbon
    if Cv < Ch
        [X C]=carve(X,C,My,'vert');
        vertRibbonCount=vertRibbonCount+1;
        %PATH(i)=2; i=i+1;
    
    %remove a horizontal ribbon
    elseif Ch < Cv
        [X C]=carve(X,C,Mx,'hor');
        horRibbonCount=horRibbonCount+1;
        %PATH(i)=1; i=i+1;
            
    %when the cost of removing a vertical and horizontal ribbon are equal
    else
        if rand < .5
            [X C]=carve(X,C,My,'vert');
            vertRibbonCount=vertRibbonCount+1;
            %PATH(i)=2; i=i+1;
        else
            [X C]=carve(X,C,Mx,'hor');
            horRibbonCount=horRibbonCount+1;
            %PATH(i)=1; i=i+1;
        end
    end
    My=cummincost(C,flex,'vert');      %compute the cumulative minimum cost
    [Cv b]=min(My(H,:,1));                 %the cost of removing a vertical ribbon
    Mx=cummincost(C,flex,'hor');
    [Ch a]=min(Mx(W,:,1));                 %the cost of removing a horizontal ribbon
    %Ch
    %a
    %Cv
    %b
end
