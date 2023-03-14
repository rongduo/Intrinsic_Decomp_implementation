function [X ] = swap( X,i,j )
%SWAP Summary of this function goes here
%   Detailed explanation goes here
    temp=X(i,:);
    X(i,:)=X(j,:);
    X(j,:)=temp;

end

