function [true_heatmap] = cal_heatmap_groundtruth(person1_data,person2_data,calib_const)
% % % Input
% person1_data: the location of the first person over time
% person2_data: the location of the second person over time
% % % Output
% true_heatmap: the true heatmap

cc = calib_const;
nX = floor((cc.x_range(2)-cc.x_range(1))/cc.x_step);
nY = floor((cc.y_range(2)-cc.y_range(1))/cc.y_step);

nSample = size(person1_data,1);
true_heatmap = zeros(nSample,nX*nY);
for i=1:nSample
    if nnz(person1_data(i,:))~=0
        x = round((person1_data(i,1)-cc.x_range(1))/cc.x_step);
        y = round((person1_data(i,2)-cc.y_range(1))/cc.y_step);
        true_heatmap(i,y*nX+x)=1;
    end
    if nnz(person2_data(i,:))~=0
        x = round((person2_data(i,1)-cc.x_range(1))/cc.x_step);
        y = round((person2_data(i,2)-cc.y_range(1))/cc.y_step);
        true_heatmap(i,y*nX+x)=1;
    end
end
end

