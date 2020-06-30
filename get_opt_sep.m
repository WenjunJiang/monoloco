function [antenna_sep,result] = get_opt_sep(csi,antenna_type,init_sep,calib_const)
%,tx_loc,rx_loc,tx_orient, rx_orient, init_sep,max_offset)

% The inital seperation measured may not be accurate. 
% The goal of this function is to adjust the separation of the antennas 
% so that the AoA or AoD calculated fit the data
%
% csi: complex csi data collected in empty space with dimension nSample*nTx*nRx*nChannel
% antenna_type: string value, it should be "tAP" if we calculate AoD, "rAP"+N if we calculate AoA
% calib_const.tx_loc: 2D coordinate of the tx with dimension nTx*2
% calib_const.rx_loc: 2D coordinate of the rx with dimension nRx*2
% calib_const.tx_orient: orientation of tx array, from tx1 to tx3, with unit rad
% calib_const.rx_orient: orientation of rx array, from rx1 to rx3, with unit rad
% init_sep: initial separation of the antennas with unit to be meter
% calib_const.max_offset: scalar, the offset we add should be in the range [0,max_offset]
% calib_const.C: speed of light
%
% antenna_sep: vector, the final seperation of the antennas
% result: output all testing results
[nChannel,nRx,nTx,nSample] = size(csi);
assert(nTx==3 && nRx==3 && nChannel==30);
cc = calib_const; % for simplicity
rAPs = {'rAP1','rAP2','rAP3'};
tAPs = {'tAP1','tAP2','tAP3'}; % Although we have only one Tx, we need an Rx to caliberate it
result = [];
if any(strcmp(rAPs,antenna_type)) % calculate AoA
    APid = str2num(antenna_type(4));
    antenna_sep = init_sep;
    % The rule is: relativeAoA = trueAoA - (rx_orient - pi)
    aoa_target = atan2(cc.tx_loc(2)-cc.rx_loc(APid,2), cc.tx_loc(1)-cc.rx_loc(APid,1))...
                    - (cc.rx_orient(APid)-pi); % the range of atan2 is [-pi,pi]
    
                
    rx2_1 = angle(csi(:,2,:,:)./csi(:,1,:,:)); % phase different between rx3 and rx1, (-pi,pi]
    rx2_1 = squeeze(rx2_1);
    d2_1 = antenna_sep(2)-antenna_sep(1);
    
    rx3_2 = angle(csi(:,3,:,:)./csi(:,2,:,:)); % phase different between rx3 and rx2, (-pi,pi]
    rx3_2 = squeeze(rx3_2);
    d3_2 = antenna_sep(3)-antenna_sep(2);
    
    rx3_1 = angle(csi(:,3,:,:)./csi(:,1,:,:)); % phase different between rx3 and rx1, (-pi,pi]
    rx3_1 = squeeze(rx3_1);
    d3_1 = antenna_sep(3)-antenna_sep(1);
    
    mse_best = 1000000;
    for offset2 = 0.001:0.001:0.05%-cc.max_offset:0.0005:cc.max_offset
%         for offset3 = -cc.max_offset:0.0005:cc.max_offset
            offset3 = 0;
            rx2_1_gt = mod(-2*pi*(d2_1+offset2)*cos(aoa_target)*cc.f(30)/cc.C+pi,2*pi)-pi;
            rx2_1_gt = repmat(rx2_1_gt',1,nChannel,nTx,nSample);
            rx2_1_gt = squeeze(rx2_1_gt);
            
            rx3_2_gt = mod(-2*pi*(d3_2-offset2+offset3)*cos(aoa_target)*cc.f(30)/cc.C+pi,2*pi)-pi;
            rx3_2_gt = repmat(rx3_2_gt',1,nChannel,nTx,nSample);
            rx3_2_gt = squeeze(rx3_2_gt);
            
            rx3_1_gt = mod(-2*pi*(d3_1+offset2+offset3)*cos(aoa_target)*cc.f(30)/cc.C+pi,2*pi)-pi;
            rx3_1_gt = repmat(rx3_1_gt',1,nChannel,nTx,nSample);
            rx3_1_gt = squeeze(rx3_1_gt);
            
            mse = (circular_mse(rx2_1,rx2_1_gt)+circular_mse(rx3_2,rx3_2_gt)+circular_mse(rx3_1,rx3_1_gt))/3;
            result = [result;offset2,offset3,circular_mse(rx2_1,rx2_1_gt),circular_mse(rx3_2,rx3_2_gt)];
            if mse<mse_best
                mse_best = mse;
                best2 = offset2;
                best3 = offset3;
            end
%         end
    end
    antenna_sep(2)=antenna_sep(2)+best2;
    antenna_sep(3)=antenna_sep(3)+best2+best3;
    antenna_sep(4)=mse_best;

elseif any(strcmp(tAPs,antenna_type)) % calculate AoA
    APid = str2num(antenna_type(4));
    antenna_sep = init_sep;
     % The rule is: relativeAoD = trueAoD - (tx_orient - pi)
    aod_target = atan2(cc.rx_loc(APid,2)-cc.tx_loc(2), cc.rx_loc(APid,1))-cc.tx_loc(1)...
                    - (cc.tx_orient-pi); % the range of atan2 is [-pi,pi]
    tx2_1 = angle(csi(:,:,2,:)./csi(:,:,1,:)); % phase different between tx3 and tx1, (-pi,pi]
    tx2_1 = squeeze(tx2_1);
    d2_1 = antenna_sep(2)-antenna_sep(1);
    
    tx3_2 = angle(csi(:,:,3,:)./csi(:,:,2,:)); % phase different between tx3 and tx2, (-pi,pi]
    tx3_2 = squeeze(tx3_2);
    d3_2 = antenna_sep(3)-antenna_sep(2);
    
    tx3_1 = angle(csi(:,:,3,:)./csi(:,:,1,:)); % phase different between tx3 and tx1, (-pi,pi]
    tx3_1 = squeeze(tx3_1);
    d3_1 = antenna_sep(3)-antenna_sep(1);
    
    mse_best = 1000000;
    for offset2 = -cc.max_offset:0.0005:cc.max_offset
        for offset3 = -cc.max_offset:0.0005:cc.max_offset
%             offset3 = -offset2;
            tx2_1_gt = mod(-2*pi*(d2_1+offset2)*cos(aod_target)*cc.f(1)/cc.C+pi,2*pi)-pi;
            tx2_1_gt = repmat(tx2_1_gt',1,nChannel,nRx,nSample);
            tx2_1_gt = squeeze(tx2_1_gt);
            
            tx3_2_gt = mod(-2*pi*(d3_2-offset2+offset3)*cos(aod_target)*cc.f(1)/cc.C+pi,2*pi)-pi;
            tx3_2_gt = repmat(tx3_2_gt',1,nChannel,nRx,nSample);
            tx3_2_gt = squeeze(tx3_2_gt);
            
            tx3_1_gt = mod(-2*pi*(d3_1+offset2+offset3)*cos(aod_target)*cc.f(1)/cc.C+pi,2*pi)-pi;
            tx3_1_gt = repmat(tx3_1_gt',1,nChannel,nRx,nSample);
            tx3_1_gt = squeeze(tx3_1_gt);
            
            mse = (circular_mse(tx2_1,tx2_1_gt)+circular_mse(tx3_2,tx3_2_gt)+circular_mse(tx3_1,tx3_1_gt))/3;
            result = [result;offset2,offset3,mse];
            if mse<mse_best
                mse_best = mse;
                best2 = offset2;
                best3 = offset3;
            end
        end
    end
    antenna_sep(2)=antenna_sep(2)+best2;
    antenna_sep(3)=antenna_sep(3)+best2+best3;
    antenna_sep(4)=mse_best;
end

end

