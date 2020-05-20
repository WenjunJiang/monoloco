function sage_const = sage_generate_steering_vector(sage_const)
% SAGE_GENERATE_STEERING_VECTOR: Generate steering vectors of time, frequency and space.
    T = sage_const.T;
    F = sage_const.F;
    A = sage_const.A;
    AoD = sage_const.D;
    frequency_interval = sage_const.frequency_interval;
    time_interval = sage_const.time_interval;
    antenna_spatio = sage_const.antenna_spatio;
    aod_spatio = sage_const.aod_spatio;
    
    doppler_range = sage_const.doppler_range;
    aoa_range = sage_const.aoa_range;
    aod_range = sage_const.aod_range;
    tof_range = sage_const.tof_range;
    
    sage_const.tof_candidates = exp(-1j * 2 * pi * frequency_interval * (0:F-1).' * tof_range);
    sage_const.aoa_candidates = exp(-1j * 2 * pi * antenna_spatio * (0:A-1).' * cos(aoa_range));
    sage_const.aod_candidates = exp(-1j * 2 * pi * aod_spatio * (0:AoD-1).' * cos(aod_range));
    sage_const.doppler_candidates = exp(1j * 2 * pi * time_interval * (0:T-1).' * doppler_range);
end