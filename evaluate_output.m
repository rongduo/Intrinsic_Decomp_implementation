clear all;

index=1:12; % List of indexes for the images that will be evaluated
% initialize a series of placeholder variables for accumulating errors from
% all of the images
sum_lmse_albedo=0;
sum_lmse_shading=0;
sum_mse_albedo=0;
sum_mse_shading=0;

% Specify the folders for ground-truth and results
address_gt='data/MPI/'; % the folder for the ground truth
address='Output_results/'; % folder for the output from the algorithm

% A placeholder array for the numerical results.
value=zeros(5*(size(index,2)+1),1);
%%
for i=index
    
    % read ground truth albedo
    a2=imread(strcat(address_gt,sprintf('gt_albedo%04d.png',i)));
    [h, w, ~]=size(a2);
    a2=im2double(a2);
    % read ground truth shading
    s2=imread(strcat(address_gt,sprintf('gt_shading%04d.png',i)));
    s2=im2double(s2);

    
    % read the albedo result
    a=imread(strcat(address,sprintf('albedo%04d.png',i))); 
    a=im2double(a);
    % read the shading result
    s=imread(strcat(address,sprintf('shading%04d.png',i)));
    s=im2double(s);
    s=s(:,:,[1 1 1]);
    
    % calculate ssim values, then transform them to DSSIM and write into
    % the csv file.
    error_mse_albedo=evaluate_one_k(a,a2);
    error_mse_shading=evaluate_one_k(s,s2);
    value(5*(i-1)+1)=roundn(error_mse_albedo,-4);
    value(5*(i-1)+2)=roundn(error_mse_shading,-4);  
    
    % calculate lmse values and write them into the xls or csv file.
    error_lmse_albedo=levaluate_one_k(a,a2);
    value(5*(i-1)+3)=roundn(error_lmse_albedo,-4);
    error_lmse_shading=levaluate_one_k(s,s2);
    value(5*(i-1)+4)=roundn(error_lmse_shading,-4);
    
       
    % accumulate the errors
    sum_mse_albedo=sum_mse_albedo+error_mse_albedo;
    sum_mse_shading=sum_mse_shading+error_mse_shading;
    sum_lmse_albedo=sum_lmse_albedo+error_lmse_albedo;
    sum_lmse_shading=sum_lmse_shading+error_lmse_shading;

    % Display the error values for current example
    disp(error_mse_albedo);
    disp(error_mse_shading);
    disp(error_lmse_albedo);
    disp(error_lmse_shading);
end
%%

%% Calculate average error values, write them into the csv file and display here
disp('Average values');
disp('MSE for albedo');
value(5*(i-0)+1)=roundn(sum_mse_albedo/size(index,2),-4);
disp(sum_mse_albedo / size(index,2));
disp('MSE for shading');
value(5*(i-0)+2)=roundn(sum_mse_shading/size(index,2),-4);
disp(sum_mse_shading/size(index,2));
disp('Average MSE for albedo and shading')
disp((sum_mse_albedo + sum_mse_shading) / 2.0 / size(index,2));


disp('LMSE for albedo');
value(5*(i-0)+3)=roundn(sum_lmse_albedo/size(index,2),-4);
disp(sum_lmse_albedo/size(index,2));
disp('LMSE for albedo');
value(5*(i-0)+4)=roundn(sum_lmse_shading/size(index,2),-4);
disp(sum_lmse_shading/size(index,2));
disp('Average LMSE for albedo and shading')
disp((sum_lmse_albedo + sum_lmse_shading) / 2.0 / size(index,2));
xlswrite([address sprintf('evaluation_results.csv')],value);