%extract the pixels specified in the mask
function imV=extract(im,mask)
	im=reshape(im,[],size(im,3));
	imV=im(find(mask),:);
end
