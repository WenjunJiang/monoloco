clc,clear;
data_path = '../../data/TRAP_050220_xy/';
file_name = 'spectrum_rAP2_jwj_hori1.mat';
load(strcat(data_path,file_name));

t=10;
figure();
subplot(2,1,1);
feature_map = reshape(tad_spectrum(t*10,:),35,25);
image(feature_map','CDataMapping','scaled');
subplot(2,1,2);
true_map = reshape(true_heatmap(t*10,:),35,25);
image(true_map','CDataMapping','scaled');
