function weighted_map= weighted_bialateral( n,edge,image,f )

[row column z]=size(image);

sigma_c=0.3;
sigma2=sigma_c*1.732*11; %color weighted
sigma3=10; %gausian weighted 
sigma4=0.2;
bsize2=2*f+1;
[A,B] = meshgrid(-f:f,-f:f);  
D = exp(-(A.^2+B.^2)/(2*sigma3^2));
sigma=zeros(bsize2*n,bsize2);
sigma=sigma(:,:,[1 1 1]);
colorkernel=zeros(bsize2*n,bsize2);
weighted_map=zeros(size(edge,1),1);
image_f=padarray(image, [f,f], 'symmetric');
position_map=zeros(n,2);
       for j=1:n
            p_x=fix(j/row)+1;
            p_y=mod(j,row);
            if(p_y==0)
                p_y = row;
                p_x=fix(j/row);
            end
            position_map(j,:)=[p_y p_x];
            sigma_p=image_f(p_y:p_y+2*f,p_x:p_x+2*f,:);
            sigma(bsize2*j-2*f:bsize2*j,:,:)=sigma_p;
            dR = abs(sigma_p(:,:,1)-image(p_y,p_x,1));
            dG = abs(sigma_p(:,:,2)-image(p_y,p_x,2));
            dB = abs(sigma_p(:,:,3)-image(p_y,p_x,3));
            color_kernel = exp(-(dR.^2+dG.^2+dB.^2)/(2*3*sigma4^2));
            colorkernel(bsize2*j-2*f:bsize2*j,:)=color_kernel.*D;
       end
       for i=1:size(edge,1)
            diff=sum(sigma(bsize2*edge(i,1)-2*f:bsize2*edge(i,1),:,:)-sigma(bsize2*edge(i,2)-2*f:bsize2*edge(i,2),:,:),3);
            kpq=exp(-sum(sum (colorkernel(bsize2*edge(i,1)-2*f:bsize2*edge(i,1),:).*(diff.^2)))/(2*3*sigma2^2));  
            weighted_map(i)=kpq;
       end
end
            






