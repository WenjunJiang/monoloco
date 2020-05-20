function [outputArg1,outputArg2] = monoloco_main(csi_data_static,csi_data_motion,cal_range)
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
% estimate, N is # of iteration, fi is interval of subcarriers, FI,TI,AS,DS
% are interval/spatio of frequency,time,AoA,AoD
% TR,AR,DR are range of ToF,AoA,Doppler
T = 100; F = 30; A = 3; AoD=3; L = 5; N = 100; FI = 2*312.5e3; TI = 1/sample_rate;
AS_Rx1 = 0.16; AS_Rx2 = 0.13; AS_Rx3 = 0.16; DS=0.07; 
G=1000; TR = (-100:1:100)*1e-9; AR = (0:1:180)/180*pi; AoDR = (0:1:180)/180*pi; DR = -20:1:20; UR = 1;
sage_const = sage_set_const(T, F, A, AoD, L, N, FI, TI, AS_Rx3, DS, TR, AR, AoDR, DR, UR);
sage_const = sage_generate_steering_vector(sage_const);

csi_data=csi_data_static(cal_range,1:270);
% Data sanitization
for jj = 1:size(csi_data,1)
    for t=1:3
        csi_data(jj,t*90-89:t*90) = csi_sanitization(csi_data(jj,t*90-89:t*90));
    end
end
% Path estimation.
estimated_parameter = sage_main(csi_data, 1000, sage_const); %5*L*G

end

