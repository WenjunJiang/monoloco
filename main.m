clear;
close all;

data_path = '../../data/csi_050220/';
param_path = '../../data/param_050220/';
if ~exist(param_path,'dir')
    mkdir(param_path)
end
motion_files = {"jwj_hori1.mat","sh_verti2.mat","shjwj_vh1.mat","jwj_hori2.mat","sh_walk1.mat","shjwj_vs1.mat","jwj_rand1.mat","sh_walk2.mat","shjwj_vv1.mat","jwj_rand2.mat","sh_walk3.mat","shjwj_vw1.mat","jwj_static1.mat","sh_walk4.mat","shjwj_vw2.mat","jwj_static2.mat","shjwj_hh1.mat","shjwj_wh1.mat","jwj_verti1.mat","shjwj_hs1.mat","shjwj_wh2.mat","jwj_verti2.mat","shjwj_hv1.mat","shjwj_ws1.mat","jwj_walk1.mat","shjwj_hw1.mat","shjwj_ws2.mat","jwj_walk2.mat","shjwj_hw2.mat","shjwj_ws3.mat","jwj_walk3.mat","shjwj_rand1.mat","shjwj_ws4.mat","jwj_walk4.mat","shjwj_rand2.mat","shjwj_wv1.mat","shjwj_rand3.mat","shjwj_wv2.mat","shjwj_sh1.mat","shjwj_ww1.mat","sh_hori1.mat","shjwj_ss1.mat","shjwj_ww2.mat","sh_hori2.mat","shjwj_ss2.mat","shjwj_ww3.mat","sh_rand1.mat","shjwj_sv1.mat","shjwj_ww4.mat","sh_rand2.mat","shjwj_sw1.mat","shjwj_ww5.mat","sh_static1.mat","shjwj_sw2.mat","shjwj_ww6.mat","sh_static2.mat","shjwj_sw3.mat","sh_verti1.mat","shjwj_sw4.mat"};

% static_file_name ='rx3_nobody_empty1.mat';
% motion_file_name ='rx3_sh_walk1.mat';
static_file_name = 'nobody_empty1.mat';
motion_file_name = 'sh_walk1.mat';

vicon_path='../vicon_segment/segment/';
vicon_file_name='sh_walk1.mat';

% set parameters
rx_loc = [3.44,0; 0,-2.5; -3.44,0];
rx_orient = [-pi/2; pi; pi/2];
tx_loc = [0,-3.48];
tx_orient = pi;
x_range=[-3.5,3.5];
y_range=[-2.5,2.5];

% % for widar part
% load([data_path,static_file_name],'csi_data');
% csi_data_static=csi_data;
% load([data_path,motion_file_name],'csi_data');
% csi_data_motion=csi_data;

for fi = 1:length(motion_files)
    motion_file_name = motion_files{fi};
    vicon_file_name = motion_files{fi};
    
    all_parameter = monoloco_main(static_file_name, motion_file_name);
    % [cal_aoa,cal_range,cal_location]=widar_main(rx3_loc,tx_loc,csi_data_static,csi_data_motion,cal_range,motion_file_name);

    % save all_parameter numRx*numParam*numPath*numSample
    save(strcat(param_path,'param_',motion_file_name),'all_parameter');

    [pos,pos_feature] = monoloco_localization(all_parameter, tx_loc, tx_orient, rx_loc, rx_orient, x_range, y_range);

    % save pos_feature %numPos*numSample*nFeature
    save(strcat(param_path,'feature_',motion_file_name),'pos','pos_feature');

%     % for vicon part
%     load(strcat(vicon_path,vicon_file_name));
%     len=size(video_data,1);
%     video_data=video_data/1000;
%     person1_data=video_data(:,19:1:21);
%     person2_data=video_data(:,70:1:72);
%     
%     pos_gt = [person1_data,person2_data];
%     save(strcat(param_path, 'pos_',vicon_file_name),'pos_gt');
%     
%     [true_aoa,true_aod, true_range]=cal_para_groundtruth(tx_loc, rx_loc,person1_data);
%     save(strcat(param_path,'p1_',vicon_file_name),'true_aoa','true_aod','true_range');
%     
%     [true_aoa,true_aod, true_range]=cal_para_groundtruth(tx_loc, rx_loc,person2_data);
%     save(strcat(param_path,'p2_',vicon_file_name),'true_aoa','true_aod','true_range');
end


% % for plot
% cal_location=cal_location.';
% figure();
% xlim(x_range);
% ylim(y_range);
% plot(person1_data(:,1), person1_data(:,2), 'r');hold on;
% % plot(person2_data(:,1), person2_data(:,2), 'b');hold on;
% plot(cal_location(:,1),cal_location(:,2),'b*-');