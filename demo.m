clear all;
tic;
run('./vlfeat-0.9.20-master/toolbox/vl_setup.m'); % setup the vlfeat library, which will be used for finding k nearest neighbors.

folder='data/MPI/';  % data_folder path
ratio=1; %resize the image to its size times ratio
outputDir = 'Output_results/';

% if the output folder does not exist, create one
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

index = 1:12; % list for the indexes of files that will be test
for id=index
	%read data
	[image depth mask chroma normals luma]=readData(folder,ratio,id);  % read data from the data folder according to index
	luma=min(max(luma(find(mask))*0.9999+0.0001,0),1);%to avoid log 0   
	points=getVectors(size(depth,1),size(depth,2),57);%horizontal field of view size
	points=points.*depth(:,:,[1 1 1]); 
    points_whitened=reshape(whiten(reshape(points,[],3)),size(image));%whiten the 3D points  
	n=numel(find(mask));%number of pixels
	log_image=max(extract(toLog(image),mask),-6);
	normals=reshape(normals,[],3);
    [row column z]=size(image);


    num_x = 10;  % parameter used for control the number of neighbors for regularizers related to albedo
    num_y = 5; % parameter used for control the number of neighbors for regularizerse related to shading
    W=sparse(1:n,1:n,sqrt(luma+0.1),n,n);
     
    
    % albedo regularization
    num_nn_albedo=num_x;
	ind=randi(n,num_nn_albedo,n);
	edge=[reshape(repmat(1:n,[num_nn_albedo 1]),[],1) reshape(ind,[],1)];%a list of edges N_A  % create random connections accross the whole image
	edge=edge(find(edge(:,1)~=edge(:,2)),:); % remove the self-connected edges
    edge=edge(find(edge(:,2)~=0),:);   
	edge=unique([min(edge(:,1),edge(:,2)) max(edge(:,1),edge(:,2))],'rows');  % If there are two edges with two same points but different orders, just keep one of them
    weight_albedo_reg = weighted_bialateral(n, edge, image,2);  % calculate the weight used for each edge within the regularizer for albedo
	num_albedo_reg=size(edge,1);
	Q1=sparse(1:num_albedo_reg,edge(:,1),weight_albedo_reg,num_albedo_reg,n)-sparse(1:num_albedo_reg,edge(:,2),weight_albedo_reg,num_albedo_reg,n);

    
    %directIrradiance regularization
	feature_directIrradiance=[reshape(points_whitened,[],3) normals]';%feature vector (x,y,z,n_x,n_y,n_z)
	num_nn_directIrradiance=num_y*2;%round(ratio*2);
	ind=double(vl_kdtreequery(vl_kdtreebuild(feature_directIrradiance),feature_directIrradiance,feature_directIrradiance,'numneighbors',num_nn_directIrradiance));
	edge=[reshape(repmat(1:n,[num_nn_directIrradiance 1]),[],1) reshape(ind,[],1)];%a list of edges N_D
	num_directIrradiance_reg=size(edge,1);
	weight_directIrradiance_reg=ones(1,num_directIrradiance_reg);  % use the same weight (1) for each pair of neighbors within the regularizer
	Q2=sparse(1:num_directIrradiance_reg,edge(:,1),weight_directIrradiance_reg,num_directIrradiance_reg,n)-sparse(1:num_directIrradiance_reg,edge(:,2),weight_directIrradiance_reg,num_directIrradiance_reg,n);
	
	%indirect Irradiance regularization
	%indirectIrradiance map should be spatially coherent
	feature_indirectIrradiance=[reshape(points,[],3)]';
	num_nn_indirectIrradiance=num_y*2;%round(ratio*2);
	ind_indirectIrradiance=double(vl_kdtreequery(vl_kdtreebuild(feature_indirectIrradiance),feature_indirectIrradiance,feature_indirectIrradiance,'numneighbors',num_nn_indirectIrradiance));
	edge=[reshape(repmat(1:n,[num_nn_indirectIrradiance 1]),[],1) reshape(ind_indirectIrradiance,[],1)];%a list of edges N_N
	num_indirectIrradiance_reg=size(edge,1);
	weight_indirectIrradiance_reg=ones(1,num_indirectIrradiance_reg); 
    Q3=sparse(1:num_indirectIrradiance_reg,edge(:,1),weight_indirectIrradiance_reg,num_indirectIrradiance_reg,n)-sparse(1:num_indirectIrradiance_reg,edge(:,2),weight_indirectIrradiance_reg,num_indirectIrradiance_reg,n);


	%indirectIrradiance map should be cloffsete to 0 in log domain
    Q4=sparse(1:n,1:n,ones(1,n),n,n);
	
	%illuminationColor regularization
	num_nn_illuminationColor=num_x;
	ind=randi(n,num_nn_illuminationColor,n);
	edge=[reshape(repmat(1:n,[num_nn_illuminationColor 1]),[],1) reshape(ind,[],1)];%a list of edges N_C
	num_illuminationColor_reg=size(edge,1);
	points_vector=reshape(points,[],3);
	diff=sum(abs(points_vector(edge(:,1),:)-points_vector(edge(:,2),:)),2);
	weight_illuminationColor_reg=(1-diff/max(diff));
    weight_illuminationColor_reg=sqrt(weight_illuminationColor_reg);
    Q5=sparse(1:num_illuminationColor_reg,edge(:,1),weight_illuminationColor_reg,num_illuminationColor_reg,n)-sparse(1:num_illuminationColor_reg,edge(:,2),weight_illuminationColor_reg,num_illuminationColor_reg,n);   
 
    %%% Start to solve the unknowns
	t=[1;0.001;10;0.1;10;5];  % List of lambdas
    log_i=reshape(log_image,[],3);
    % initialize each components
    reflectance=zeros(n,3);
    direct=zeros(n,1);
    indirect=zeros(n,1);
    color=zeros(n,3);

    
    s=size(image);
    [reflectance  direct indirect color ]= solver(n,s,log_i,t(2),t(3),t(4),t(5),t(6),t(1),W, Q1,Q2,Q3,Q4,Q5,reflectance,direct,indirect,color); % solve each component in an iterative manner
    % transform the components from log domain to the normal range
    reflectance=exp(reshape(reflectance,size(image)));
    direct=exp(direct);
    direct=reshape(repmat(direct,1,3),size(image));
    indirect=exp(reshape(indirect,size(image,1),size(image,2)));
    indirect=repmat(indirect,[1 1 3]);
    color=exp(reshape(color,size(image)));
    % calculate the final shading
    shading=direct.*indirect.*color;
    [reflectance shading]=adjust(reflectance, shading);
    direct=shading./(color.*indirect);

	%save the output image
	imwrite(reflectance,[outputDir sprintf('albedo%04d.png',id)]);
	imwrite(shading,[outputDir sprintf('shading%04d.png',id)]);			
end
toc;