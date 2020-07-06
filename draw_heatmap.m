clc,clear;
data_path = '../../data/TRAP_050220_xy/';
file_name = 'jwj_hori1.mat';
dt = 0.1; % set delay time to 0.1s
file_AP1 = strcat('spectrum_rAP1_',file_name);
file_AP2 = strcat('spectrum_rAP2_',file_name);
file_AP3 = strcat('spectrum_rAP3_',file_name);

data1 = load(strcat(data_path,file_AP1));
data2 = load(strcat(data_path,file_AP2));
data3 = load(strcat(data_path,file_AP3));

%draw gif
h=figure;
axis tight manual % this ensures that getframe() returns a consistent size
[file_path,name,ext] = fileparts(file_name);
file_gif = strcat(name,'.gif');

for t=1:60
subplot(2,2,1);
true_map = reshape(data1.true_heatmap(t*10,:),35,25);
image(true_map','CDataMapping','scaled');

subplot(2,2,2);
feature_map = reshape(data1.tad_spectrum(t*10,:),35,25);
image(feature_map','CDataMapping','scaled');

subplot(2,2,3);
feature_map = reshape(data2.tad_spectrum(t*10,:),35,25);
image(feature_map','CDataMapping','scaled');

subplot(2,2,4);
feature_map = reshape(data3.tad_spectrum(t*10,:),35,25);
image(feature_map','CDataMapping','scaled');
drawnow

% Capture the plot as an image 
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

% Write to the GIF File 
if t == 1 
  imwrite(imind,cm,file_gif,'gif', 'writeMode','overwrite','Loopcount',inf,'delaytime',dt); 
else 
  imwrite(imind,cm,file_gif,'gif','WriteMode','append','delaytime',dt); 
end

end