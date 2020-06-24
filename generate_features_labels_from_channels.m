%% Script that converts channels to features and labels for Dloc
% loads the dataset from DATASET_NAME and saves the datasets
clearvars
%% Tuneable Parameters
DATASET_NAME = 'July16'; %'July16','July18','July18_different_APsHeight','July22_1_ref','July22_2_ref','jacobs_July28','jacobs_July28_2'
ANT_SEP = 0.0259; % antenna separation on the access points
OUTPUT_GRID_SIZE = 0.1; %the output grid size of each pixel
OUTPUT_SIGMA = 0.25; % the gaussian variance of the ouput gaussian target
%%
load(['channels_,',DATASET_NAME,'.mat']);
[n_points,n_sub,n_ant,n_ap] = size(channels);

%% Estimating AoA
S = get_2dsteering_matrix(theta_vals,d_vals,n_sub,mean(opt.freq),...
    mean(diff(opt.freq)), opt.ant_sep);
parfor i=1:n_points
    [aoa_pred(i,:),d_pred(i,:)] = get_least_tofs_aoa(...
        squeeze(channels(i,:,:,:)),theta_vals,d_vals,opt,n_ap,S);
    if(mod(i,1000)==0)
        disp(i)
    end
end
%% Real ToF compensation
for i=1:n_points
    parfor j=1:n_ap
        channels_wo_offset(i,:,:,j) = squeeze(channels(i,:,:,j)).*...
            exp( 1j*2*pi*opt.freq.'*( d_pred(i,j) - real_tof(i,j) )./3e8 );
    end
end
%% 

max_x = 3*d1(end)/2;
max_y = 3*d2(end)/2;
min_x = -d1(end)/2;
min_y = -d2(end)/2;
d1 = min_x:input_grid_size:max_x;
d2 = min_y:input_grid_size:max_y;

features_w_offset = zeros(n_points,n_ap,length(d2),length(d1));
features_wo_offset = zeros(n_points,n_ap,length(d2),length(d1));

%% Get features
parfor i=1:n_points
    features_with_offset(i,:,:,:) = generate_features_abs(squeeze(...
        channels(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
    features_without_offset(i,:,:,:) = generate_features_abs(squeeze(...
        channels_wo_offset(i,:,:,:)),ap,theta_vals,d_vals,d1,d2,opt);
    if(mod(i,1000)==0)
        disp(i);
    end    
end
labels_gaussian_2d = get_gaussian_labels(labels,...
    OUTPUT_GRID_SIZE,output_sigma,d1,d2);
save(['dataset_',DATASET_NAME,'.mat'], 'features_with_offset',...
    'features_without_offset','labels_gaussian_2d','labels','-v7.3');