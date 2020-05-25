function [all_parameter] = monoloco_main(static_file_name,motion_file_name)
% set parameters
data_path = '../../data/csi_050220/';
save_path = '../Temp/';
c=299792458;
carrier_frequency=5825000000.00000;
sample_rate=1000;

xb=[-3.5,3.5];
yb=[-3.5,3.5];

% SAGE parameters
% T is length of segments, F is # of subcarrier, A is # of receiving antennas, L is # of paths to
% estimate, N is # of iteration, fi is interval of subcarriers, FI,TI,AS,DS
% are interval/spatio of frequency,time,AoA,AoD
% TR,AR,DR are range of ToF,AoA,Doppler
T = 100; F = 30; A = 3; AoD=3; L = 5; N = 100; FI = 2*312.5e3; TI = 1/sample_rate;
AS = [0.16, 0.13, 0.16]; DS=0.07;
G=1000; TR = (-100:1:100)*1e-9; AR = (0:1:180)/180*pi; AoDR = (0:1:180)/180*pi; DR = -20:1:20; UR = 1;

numRx = length(AS);
all_parameter = zeros(numRx,5,L,G); %estimated parameter of all receivers 3*5*L*G
for rxId = 1:length(AS)
    sage_const = sage_set_const(T, F, A, AoD, L, N, FI, TI, AS(rxId), DS, TR, AR, AoDR, DR, UR);
    sage_const = sage_generate_steering_vector(sage_const);

    data = load(strcat(data_path,'rx',num2str(rxId),'_',motion_file_name));
    csi_data=data.csi_data(:,1:270);
    % Data sanitization
    for jj = 1:size(csi_data,1)
        for t=1:3
            csi_data(jj,t*90-89:t*90) = csi_sanitization(csi_data(jj,t*90-89:t*90));
        end
    end
    % Path estimation.
    estimated_parameter = sage_main(csi_data, G, sage_const); %5*L*G
    
    G2 = size(estimated_parameter,3);
    if G2<size(all_parameter,4)
        all_parameter = all_parameter(:,:,:,1:G2);
    end
    all_parameter(rxId,:,:,:)=estimated_parameter;
end


end

