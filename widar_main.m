function [aoa,range,location] = widar_main(rx_loc,tx_loc,csi_data_static,csi_data_motion,cal_range,trace_file)
% main function for widar cal
% input rx_loc 

% set parameters
data_path = '../csi_segment/segment/';
save_path = '../Temp/';
c=299792458;
carrier_frequency=5825000000.00000;
sample_rate=1000;

xb=[-3.5,3.5];
yb=[-3.5,3.5];

% SAGE parameters
% T is length of segments, F is # of subcarrier, A is # of receiving antennas, L is # of paths to
% estimate, N is # of iteration, fi is interval of subcarriers, FI,TI,AS
% are interval/spatio of frequency,time,AoA
% TR,AR,DR are range of ToF,AoA,Doppler
T = 100; F = 30; A = 3; L = 5; N = 100; FI = 2*312.5e3; TI = 1/sample_rate; AS = 0.5; G=1000;
TR = (-100:1:100)*1e-9; AR = (0:1:180)/180*pi; DR = -20:1:20; UR = 1;
sage_const = sage_set_const(T, F, A, L, N, FI, TI, AS, TR, AR, DR, UR);
sage_const = sage_generate_steering_vector(sage_const);

% Path matching parameters
MN = round(sample_rate / T) + 1;
MS = 3;
parameter_weight = [1e9,180/pi,1,1].' .* [1/200,1/90,1/80,10].';

if ~exist([save_path,'device_orient.mat'], 'file')

    % % Load static data
    % spth = [data_path,static_file];
    % load(spth);

    csi_data=csi_data_static(cal_range,1:90);
    % csi_data=csi_data(1:10000,1:90);
    % csi_data=csi_data(:,1:90);

    % Data sanitization
    for jj = 1:size(csi_data,1)
        csi_data(jj,:) = csi_sanitization(csi_data(jj,:));
    end

    % Path estimation.
    estimated_parameter = sage_main(csi_data, G, sage_const);
%     save([save_path,'PPM-',static_file],'estimated_parameter');

    % Path matching
    estimated_path = sage_path_mapping_filter(estimated_parameter, parameter_weight, MN ,MS);
%     save([save_path,'PMM-',static_file],'estimated_path');

    % Identify LoS path
    [~, los_index] = max(sum(abs(estimated_path(4,:,:)),3));
    los_path = squeeze(estimated_path(:,los_index,:));
    los_aoa = median(los_path(2,:)); %\phi_Tx in Figure 10

    % Calculate device orientation
    rx_rmat = [cos(los_aoa), sin(los_aoa); -sin(los_aoa), cos(los_aoa)]; % rotate anti-clockwise by los_aoa
    rx_rvec = (rx_loc-tx_loc) / norm(rx_loc-tx_loc); % I think here "rx_loc" should be "rx_loc - tx_loc"
    rx_dvec = rx_rvec.' * rx_rmat;
    orient = atan2(rx_dvec(2), rx_dvec(1)); %\psi_r in Figure 10
        
    save([save_path,'device_orient.mat'],'orient'); % one for each receiver!
else
    load([save_path,'device_orient.mat']);
end

csi_data=csi_data_motion(cal_range,1:90);

% temp_switch_2=false;
% if temp_switch_2
if ~exist([save_path,'PPM-',trace_file], 'file')
    % Find reference antenna.
    csi_amplitude = mean(abs(csi_data));
    csi_variance = sqrt(var(abs(csi_data)));
    csi_ratio = csi_amplitude ./ csi_variance;
    ant_ratio = mean(reshape(csi_ratio, F, A), 1);
    [~, midx] = max(ant_ratio);
    csi_ref = repmat(csi_data(:,(midx-1)*F+(1:F)), 1, A);

    % Weight
    alpha_all = 0;
    for jj = 1:size(csi_data, 2)
        alpha = min(abs(csi_data(:,jj)));
        alpha_all = alpha_all + alpha;
        csi_data(:,jj) = (abs(csi_data(:,jj)) - alpha) .* exp(1j * angle(csi_data(:,jj)));
    end
    beta = alpha_all / size(csi_data,2) * 1000;
    for jj = 1:size(csi_ref,2)
        csi_ref(:,jj) = (abs(csi_ref(:,jj)) + beta) .* exp(1j * angle(csi_ref(:,jj)));
    end
    csi_mult = csi_data .* conj(csi_ref);

    % Filter
    hlfrt = sample_rate / 2;
    uppe_orde = 6;
    uppe_stop = 80;
    lowe_orde = 3;
    lowe_stop = 2;
    [B,A] = butter(5, [lowe_stop/hlfrt, uppe_stop/hlfrt], 'bandpass');
    csi_filter = zeros(size(csi_mult));
    for kk = 1:size(csi_mult,2)
        csi_filter(:,kk) = filtfilt(B,A,csi_mult(:,kk));
    end

    % Path estimation
    estimated_parameter = sage_main(csi_filter, G, sage_const);
    save([save_path,'PPM-',trace_file],'estimated_parameter');
else
    load([save_path,'PPM-',trace_file]);
end

% Path matching
if ~exist([save_path,'PMM-',trace_file], 'file')
    estimated_path = sage_path_mapping_filter(estimated_parameter, parameter_weight, MN ,MS);
    save([save_path,'PMM-',trace_file],'estimated_path');
else
    load([save_path,'PMM-',trace_file]);
end

% Identify reflection
[~, rfl_index] = max(sum(abs(estimated_path(4,:,:)),3));
rfl_path = squeeze(estimated_path(:,rfl_index,:));

% Remove first and last second noisy measurements
motion_index = 11:size(rfl_path,2)-10;

% Parameter denoise.
tof = hampel(rfl_path(1,motion_index),round(0.5*(sample_rate/T)),2e-9).';
aoa = hampel(rfl_path(2,motion_index),round(0.5*(sample_rate/T)),pi/18).';
aoa = smooth(aoa, round(0.5*(sample_rate/T)));
doppler = hampel(rfl_path(3,motion_index).', round(0.5*(sample_rate/T)),6);
doppler = smooth(doppler, round(0.5*(sample_rate/T)));
power = rfl_path(4,motion_index).';

% Kalman Filter
range = sage_tof_filter(tof*c, -doppler*c/carrier_frequency, T/sample_rate).';
range(range <= 0.3) = 0.3;
range = smooth(range, round(sample_rate/T));

% range = smooth(range, round(0.5*(sample_rate/T)));

% Transformation
range = range + norm(rx_loc - tx_loc);
aoa = orient - aoa;

% Localization
location = sage_localization(range, aoa, rx_loc, xb, yb);
for ii = 1:size(location,1)
    location(ii,:) = smooth(location(ii,:), round(0.5*(sample_rate/T)));
end

end