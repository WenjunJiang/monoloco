clear;
close all;
clc;

data_path = '../../data/csi_050220/';
param_path = '../../data/param_050220/';
vicon_path='../vicon_segment/segment/';
spectrum_path = '../../data/spectrum_050220/';
if ~exist(param_path,'dir')
    mkdir(param_path)
end
if ~exist(spectrum_path,'dir')
    mkdir(spectrum_path)
end
motion_files = {"sh_verti2.mat","shjwj_vh1.mat","jwj_hori2.mat","sh_walk1.mat","shjwj_vs1.mat","jwj_rand1.mat","sh_walk2.mat","shjwj_vv1.mat","jwj_rand2.mat","sh_walk3.mat","shjwj_vw1.mat","jwj_static1.mat","sh_walk4.mat","shjwj_vw2.mat","jwj_static2.mat","shjwj_hh1.mat","shjwj_wh1.mat","jwj_verti1.mat","shjwj_hs1.mat","shjwj_wh2.mat","jwj_verti2.mat","shjwj_hv1.mat","shjwj_ws1.mat","jwj_walk1.mat","shjwj_hw1.mat","shjwj_ws2.mat","jwj_walk2.mat","shjwj_hw2.mat","shjwj_ws3.mat","jwj_walk3.mat","shjwj_rand1.mat","shjwj_ws4.mat","jwj_walk4.mat","shjwj_rand2.mat","shjwj_wv1.mat","shjwj_rand3.mat","shjwj_wv2.mat","shjwj_sh1.mat","shjwj_ww1.mat","sh_hori1.mat","shjwj_ss1.mat","shjwj_ww2.mat","sh_hori2.mat","shjwj_ss2.mat","shjwj_ww3.mat","sh_rand1.mat","shjwj_sv1.mat","shjwj_ww4.mat","sh_rand2.mat","shjwj_sw1.mat","shjwj_ww5.mat","sh_static1.mat","shjwj_sw2.mat","shjwj_ww6.mat","sh_static2.mat","shjwj_sw3.mat","sh_verti1.mat","shjwj_sw4.mat"};
rx_loc = [3.44,0; 0,-2.5; -3.44,0];
rx_orient = [-pi/2; pi; pi/2];
tx_loc = [0,-3.48];
tx_orient = pi;
x_range=[-3.5,3.5];
y_range=[-2.5,2.5];

c=299792458;
% c = 3e8;
carrier_frequency=5825000000.00000;
sample_rate=1000;

T = 100; F = 30; A = 3; AoD=3; L = 5; N = 100; FI = 2*312.5e3; TI = 1/sample_rate;
AS_vec = [0,0.016,0.032;...
          0,0.013,0.026;...
          0,0.016,0.032]; % 1D distance, num of rAP * num of antennas on each AP
DS_vec=[0,0.0085,0.0105];
G=10; ToFR = (1:1:100)*1e-9; AoAR = (2:2:180)/180*pi; AoDR = (2:2:180)/180*pi; DR = -20:1:20; UR = 1;
ToFR=0;
max_offset = 0.01; %The offset we can compensate antenna

needCalib = true;

calib_const = calib_set_const(tx_loc,rx_loc,tx_orient,rx_orient,c,carrier_frequency,FI,max_offset);
tic;

if needCalib
    for rAPId = [1,3]
        data_empty = load(strcat(data_path,'rx',num2str(rAPId),'_nobody_empty1.mat'));
        csi_data_empty=data_empty.csi_data(:,1:270);
        %% Data sanitization
%         for jj = 1:size(csi_data_empty,1)
%             for t=1:3
%                 csi_data_empty(jj,t*90-89:t*90) = csi_sanitization(csi_data_empty(jj,t*90-89:t*90));
%             end
%         end
        csi_data_empty = reshape(csi_data_empty',F,A,AoD,[]); 
%         csi_data_empty = permute(csi_data_empty,[4,3,2,1]); % reshape the data to nSample*nTx*nRx*nChannel
        
        [antenna_sep,result] = get_opt_sep(csi_data_empty,strcat('rAP',num2str(rAPId)),AS_vec(rAPId,:),calib_const);
        disp(strcat('rAP',num2str(rAPId),':',num2str(antenna_sep)));
        disp(result);
%         [antenna_sep_tx,result] = get_opt_sep(csi_data_empty,strcat('tAP',num2str(rAPId)),DS_vec,calib_const);
%         disp(strcat('tAP',num2str(rAPId),':',num2str(antenna_sep_tx)));
%         disp(result);
    end
end

for fi = 1:length(motion_files)

    motion_file_name = motion_files{fi};
    vicon_file_name = motion_files{fi};
    
    %     % for vicon part
    load(strcat(vicon_path,vicon_file_name));
    len=size(video_data,1);
    video_data=video_data/1000;
    person1_data=video_data(:,19:1:21);
    person2_data=video_data(:,70:1:72);
    
    [true_aoa1,true_aod1, true_range1]=cal_para_groundtruth(tx_loc, rx_loc,person1_data);
    [true_aoa2,true_aod2, true_range2]=cal_para_groundtruth(tx_loc, rx_loc,person2_data);
    
    true_aoa1_relative = rad2deg(true_aoa1-repmat(rx_orient'-pi,size(true_aoa1,1),1)); % measure subject1
    true_aoa2_relative = rad2deg(true_aoa2-repmat(rx_orient'-pi,size(true_aoa2,1),1)); % measure subject2
    aoa_tx_relative = rad2deg(atan2(tx_loc(2)-rx_loc(:,2),tx_loc(1)-rx_loc(:,1))-(rx_orient-pi))'; % measure Tx
    
    numRAP = size(AS_vec,1);
    all_parameter = zeros(numRAP,5,L,G); %estimated parameter of all receivers 3*5*L*G
    
%     for rAPId = 1:numRAP
%         sage_const = sage_set_const(T, F, A, AoD, L, N, FI, TI, AS_vec(rAPId,:), DS_vec, ToFR, AoAR, AoDR, DR, UR,c);
%         sage_const = sage_generate_steering_vector(sage_const); 

%         data = load(strcat(data_path,'rx',num2str(rAPId),'_',motion_file_name));
%         csi_data=data.csi_data(:,1:270);
%         %% Data sanitization
%         for jj = 1:size(csi_data,1)
%             for t=1:3
%                 csi_data(jj,t*90-89:t*90) = csi_sanitization(csi_data(jj,t*90-89:t*90));
%             end
%         end
%         %% Estimate 3 profiles AoA-AoD, AoA-ToF, AoD-ToF%
%         aTot = gridVecBackscatter3D(ToFR,AoAR,AoDR,AS_vec(rAPId,:),DS_vec,calib_const,true); % Get steering vector of each possible combination
%         lenSeg = 100; %We calculate one spectrum every 0.1s
%         tad_spectrum = []; % nSAmple_nChannel*nRx*nTx
%         for i=1:lenSeg:size(csi_data,1)-lenSeg+1
%             disp(i);
%             toc;
%     %         csi_data1 = reshape(csi_data(i:i+lenSeg,:)',F,A,AoD,[]); % reshape to nChannel*nRx*nTx*nSample
%     %         csi_data1 = permute(csi_data1,[2,3,1,4]); % permute to nRx*nTx*nChannel*nSample
%     %         csi_data1 = reshape(csi_data1,AoD*A,[]); % simplified version of "formatCSI"
%             csi_data1 = csi_data(i:i+lenSeg-1,:)';% simplified version of "formatCSI",nChannel*nRx*nTx_nSample
% 
%             nComps = 2; % 2 paths
%             [Pn,Ps,Qn,Qs,EigenInfo] = GetQnBackscatter(csi_data1, nComps);
%             [delayLoS,aoaLoS,aodLoS]= GetLoSParam(strcat('rAP',num2str(rAPId)),calib_const);
%             useNoise = 0;
%             music_spectrum = RAPMusicGrid(aTot,Qn,Qs,delayLoS,aoaLoS,aodLoS,AS_vec(rAPId,:),DS_vec,useNoise,calib_const);
%             tad_spectrum = [tad_spectrum;music_spectrum];
%         end
%         save(strcat(spectrum_path,'spectrum_rAP',num2str(rAPId),'_',motion_file_name),'tad_spectrum', ...
%             'ToFR','AoAR','AoDR');
%     end
    %%%%%%%%%%%%%%%%%%%%%%
%     csi_data = reshape(csi_data',F,A,AoD,[]); % reshape to nChannel*nRx*nTx*nSample
%     csi_data = permute(csi_data,[4,3,2,1]); % reshape the data to nSample*nTx*nRx*nChannel
%     
%     %%% caliberate aoa 
% %     %%%%% estimation 5 paths using SAGE %%%%%
% %     estimated_parameter = sage_main(csi_data, G, sage_const); %5*L*G
% % 
% %     G2 = size(estimated_parameter,3);
% %     if G2<size(all_parameter,4)
% %         all_parameter = all_parameter(:,:,:,1:G2);
% %     end
% %     all_parameter(rxId,:,:,:)=estimated_parameter;
%     %%%% estimate 1 path using MUSIC %%%%
%     %% Estimating AoA
%     n_points = G;
%     S = get_2dsteering_matrix(sage_const.aoa_range,sage_const.tof_range,sage_const.F,...
%         carrier_frequency,FI, AS_vec(rxId,2:end),sage_const.A,sage_const.c);
%     aoa_pred = zeros(n_points,sage_const.D); %Here we treat different Tx antenna as different AP
%     d_pred = zeros(n_points,sage_const.D);
%     for i=1:n_points
%         sample = reshape(csi_data(i,:)',sage_const.F,sage_const.A,sage_const.D); 
%         sample = permute(sample,[3,2,1]);%nTx * nRx * nChannel
%         assert(all(size(sample)==[sage_const.D,sage_const.A,sage_const.F]));
%         AoASample = permute(sample,[3,2,1]); % Here we treat different Tx antenna as different AP
%         [aoa_pred(i,:),d_pred(i,:)] = get_least_tofs_aoa(...
%             AoASample,sage_const.aoa_range,sage_const.tof_range,1.0,sage_const.D,S);
%         if(mod(i,1)==0)
%             disp(i)
%         end
%     end

end