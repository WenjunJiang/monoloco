function estimated_parameter = sage_main(csi_data, G, sage_const)
% SAGE_MAIN: Segment CSI data for parameter estimation.
    T = sage_const.T;
    L = sage_const.L;
    F = sage_const.F;
    A = sage_const.A;
    AoD = sage_const.D;
%     G = min(G, floor(size(csi_data,1)/T));
    GI = 1:T:(G-1)*T+1; % no overlap
    estimated_parameter = zeros(5, L, G);
    initial_parameter = zeros(5,L);
    initial_index = zeros(4,L);
    initial_index(1,:) = (0 - sage_const.tof_range(1)) / (sage_const.tof_range(2)-sage_const.tof_range(1)) + 1; %index of 0 in tof_range
    initial_index(2,:) = (0 - sage_const.aoa_range(1)) / (sage_const.aoa_range(2)-sage_const.aoa_range(1)) + 1;
    initial_index(3,:) = (0 - sage_const.aod_range(1)) / (sage_const.aod_range(2)-sage_const.aod_range(1)) + 1;
    initial_index(4,:) = (0 - sage_const.doppler_range(1)) / (sage_const.doppler_range(2)-sage_const.doppler_range(1)) + 1;
    residue_errors = zeros(1,G);
    fprintf('SAGE Main...');
    for ii = 1:G
        fprintf('%d,',ii);
        csi_signal = csi_data(GI(ii)+(0:T-1),:).';
        csi_signal = permute(reshape(csi_signal, F, A, AoD, T), [4,1,2,3]);
        [estimated_parameter(:,:,ii),residue_errors(:,ii)] = sage_sfg(csi_signal, initial_parameter, initial_index, sage_const);
    end
    fprintf('\n');
end