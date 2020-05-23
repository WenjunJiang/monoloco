clear;
close all;

data_path = '../../20200502/csi_050220/';

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

cal_range=1:10000;

% % for widar part
% load([data_path,static_file_name],'csi_data');
% csi_data_static=csi_data;
% load([data_path,motion_file_name],'csi_data');
% csi_data_motion=csi_data;

all_parameter = monoloco_main(static_file_name, motion_file_name,cal_range);
% [cal_aoa,cal_range,cal_location]=widar_main(rx3_loc,tx_loc,csi_data_static,csi_data_motion,cal_range,motion_file_name);

[pos,pos_prob] = monoloco_localization(all_parameter, tx_loc, tx_orient, rx_loc, rx_orient, x_range, y_range);

% for vicon part
load([vicon_path,vicon_file_name]);
len=size(video_data,1);
video_data=video_data/1000;
person1_data=video_data(:,19:1:21);
person2_data=video_data(:,70:1:72);

[true_aoa,true_range]=cal_para_groundtruth(rx3_loc,person1_data);

% for plot
cal_location=cal_location.';
figure();
xlim(x_range);
ylim(y_range);
plot(person1_data(:,1), person1_data(:,2), 'r');hold on;
% plot(person2_data(:,1), person2_data(:,2), 'b');hold on;
plot(cal_location(:,1),cal_location(:,2),'b*-');