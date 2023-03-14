%Normalized RGB Map
function [nrgb k]=NRGBMap(im)
	im=im+1e-6;
	k=sum(im,ndims(im));
    nrgb=im./repmat(k,[ones(1,ndims(im)-1) 3]);
end