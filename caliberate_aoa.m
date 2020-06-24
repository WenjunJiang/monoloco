clear;
close all;

data_path = '../../data/csi_050220/';
param_path = '../../data/param_050220/';
vicon_path='../vicon_segment/segment/';
if ~exist(param_path,'dir')
    mkdir(param_path)
end
motion_files = {"jwj_hori1.mat"};
rx_loc = [3.44,0; 0,-2.5; -3.44,0];
rx_orient = [-pi/2; pi; pi/2];
tx_loc = [0,-3.48];
tx_orient = pi;
x_range=[-3.5,3.5];
y_range=[-2.5,2.5];

c=299792458;
carrier_frequency=5825000000.00000;
sample_rate=1000;

T = 100; F = 30; A = 3; AoD=3; L = 5; N = 100; FI = 2*312.5e3; TI = 1/sample_rate;
AS = [0.16, 0.13, 0.16]; DS=0.07;
G=10; TR = (-100:1:100)*1e-9; AR = (0:1:180)/180*pi; AoDR = (0:1:180)/180*pi; DR = -20:1:20; UR = 1;

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
    
    true_aoa1_relative = rad2deg(true_aoa1-repmat(rx_orient'-pi,size(true_aoa1,1),1));
    true_aoa2_relative = rad2deg(true_aoa2-repmat(rx_orient'-pi,size(true_aoa2,1),1));
    aoa_tx_relative = rad2deg(atan2(tx_loc(2)-rx_loc(:,2),tx_loc(1)-rx_loc(:,1))-(rx_orient-pi))';
    
    numRx = length(AS);
    all_parameter = zeros(numRx,5,L,G); %estimated parameter of all receivers 3*5*L*G
    
    rxId = 1;
    sage_const = sage_set_const(T, F, A, AoD, L, N, FI, TI, AS(rxId), DS, TR, AR, AoDR, DR, UR);
    sage_const = sage_generate_steering_vector(sage_const);

    data = load(strcat(data_path,'rx',num2str(rxId),'_',motion_file_name));
    csi_data=data.csi_data(:,1:270);
    % Data sanitization
%     for jj = 1:size(csi_data,1)
%         for t=1:3
%             csi_data(jj,t*90-89:t*90) = csi_sanitization(csi_data(jj,t*90-89:t*90));
%         end
%     end
%     %%%%% estimation 5 paths using SAGE %%%%%
%     estimated_parameter = sage_main(csi_data, G, sage_const); %5*L*G
% 
%     G2 = size(estimated_parameter,3);
%     if G2<size(all_parameter,4)
%         all_parameter = all_parameter(:,:,:,1:G2);
%     end
%     all_parameter(rxId,:,:,:)=estimated_parameter;
    %%%% estimate 1 path using MUSIC %%%%
    %% Estimating AoA
    n_points = G;
    S = get_2dsteering_matrix(sage_const.aoa_range,sage_const.tof_range,sage_const.F,...
        carrier_frequency,FI, AS(rxId),sage_const.A);
    aoa_pred = zeros(n_points,sage_const.D); %Here we treat different Tx antenna as different AP
    d_pred = zeros(n_points,sage_const.D);
    for i=1:n_points
        sample = reshape(csi_data(i,:),sage_const.D,sage_const.A,sage_const.F); %nTx * nRx * nChannel
        assert(all(size(sample)==[sage_const.D,sage_const.A,sage_const.F]));
        AoASample = permute(sample,[3,2,1]); % Here we treat different Tx antenna as different AP
        [aoa_pred(i,:),d_pred(i,:)] = get_least_tofs_aoa(...
            AoASample,sage_const.aoa_range,sage_const.tof_range,1.0,sage_const.D,S);
        if(mod(i,1)==0)
            disp(i)
        end
    end

end