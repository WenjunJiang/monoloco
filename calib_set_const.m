function [calib_const] = calib_set_const(tx_loc,rx_loc,tx_orient,rx_orient,C,carrier_frequency,FI,max_offset)
% set parameter for calib
    f = (-14.5:1:14.5)*FI+carrier_frequency;
    calib_const = struct('tx_loc',tx_loc,'rx_loc',rx_loc,'tx_orient',tx_orient, ...
        'rx_orient',rx_orient,'C',C,'f',f,'fc',carrier_frequency,'max_offset',max_offset);
end

