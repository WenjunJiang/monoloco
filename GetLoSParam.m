function [delayLoS,aoaLoS,aodLoS] = GetLoSParam(antenna_type,calib_const)
% % %Input
% antenna_type: string value, it should be "rAP"+N if we calculate AoA
% % %Output
% delayLoS: delay of LoS path
% aoaLoS: angle of arrival of LoS path
% aodLoS: angle of departure of LoS path
rAPs = {'rAP1','rAP2','rAP3'};
cc = calib_const; % for simplicity
if any(strcmp(rAPs,antenna_type)) % calculate AoA
    APid = str2num(antenna_type(4));
    aoaLoS = atan2(cc.tx_loc(2)-cc.rx_loc(APid,2), cc.tx_loc(1)-cc.rx_loc(APid,1))...
                    - (cc.rx_orient(APid)-pi); % the range of atan2 is [-pi,pi]
    aodLoS = atan2(cc.rx_loc(APid,2)-cc.tx_loc(2), cc.rx_loc(APid,1))-cc.tx_loc(1)...
                    - (cc.tx_orient-pi); % the range of atan2 is [-pi,pi]
    delayLoS = sqrt(sum((cc.tx_loc-cc.rx_loc(APid,:)).^2))/cc.C; % Time of flight
end
end

