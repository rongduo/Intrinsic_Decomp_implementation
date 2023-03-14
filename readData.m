%read the data. image is the rgb image. Z is the depth map. R is the chromaticity map. N is the normal map and B is the luminance map, defined as (R+G+B)/3.
function [image Z mask R N B]=readData(address,ratio,id)
	image=im2double(imread(strcat(address,sprintf('image%04d.png',id))));
	Z=load(strcat(address,sprintf('Z%04d.mat',id)));
	
	if(exist(strcat(address,'mask.png'),'file'))
		mask=imread(strcat(address,'mask.png'))>0.5;
		mask=mask(:,:,1);	
	else
		mask=ones(size(image,1),size(image,2));  
	end
	
	Z=double(Z.Z); 
	Z(Z>10000)=-1;
	Z(Z==-1)=max(max(Z))+1;
	Z(isnan(Z))=0;          
	Z(find(mask==0))=0;
	Z=smooth(image,Z,mask);
	
	R=NRGBMap(image);%R--->chromaticity map
	
	B=mean(image,3);   %B--->luminance map  
	B=min(max(imresize(B,ratio),0),1);  
	N=getNormal(Z);  
	N=imresize(N,ratio);
	N=N./repmat(sqrt(sum(N.^2,3)),[1 1 3]);
	image=min(max(imresize(image,ratio),0),1);  
	Z=imresize(Z,ratio);
	R=imresize(R,ratio);
	mask=imresize(double(mask),ratio)>0.5;  
end
