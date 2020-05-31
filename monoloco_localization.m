function [pos, pos_feature] = monoloco_localization(all_parameter, tx_loc, tx_orient, rx_loc, rx_orient, x_range, y_range)
c=299792458;
numRx = size(all_parameter,1);
numPath = size(all_parameter,3);
G = size(all_parameter,4);
varTan1=0.2; varTan2=0.2; varDist=1;
x_step = 0.2; y_step = 0.2; % eash square is 20 cm by 20 cm
pos = []; 
for i = x_range(1):x_step:x_range(2)
    for j = y_range(1):y_step:y_range(2)
        pos = [pos;i,j];
    end
end
numPos = size(pos,1);
pos_pdf = ones(numRx,numPath-1,numPos,3,G); % numRx*(numPath-1)*numPos*numPdf*G
    
for rxId = 1:numRx
    los_param = squeeze(all_parameter(rxId,:,1,:)); % 5*G we treat the first path as LOS
    tof_tx = los_param(1,:); aoa_tx = los_param(2,:); aod_tx = los_param(3,:);
    
    rx_loc_n = rx_loc(rxId,:);
    rx_orient_n = rx_orient(rxId);
    for pId = 2:numPath
        path_param = squeeze(all_parameter(rxId,:,pId,:));
        tof_rx = path_param(1,:); aoa_rx = path_param(2,:); aod_rx = path_param(3,:);
        for posId = 1:size(pos,1)
            x = pos(posId,1);
            y = pos(posId,2);

            dist1 = abs(atan2(y-tx_loc(2), x-tx_loc(1)) - (tx_orient-pi+aod_rx));
            dist2 = abs(atan2(y-rx_loc_n(2), x-rx_loc_n(1)) - (rx_orient_n+aoa_rx-pi));
            dist3 = abs(norm([x,y]-rx_loc_n) + norm([x,y]-tx_loc) - norm(rx_loc_n-tx_loc) - (tof_rx-tof_tx)*c);            
%             dist1 = abs((y-tx_loc(2))/(x-tx_loc(1)) - tan(tx_orient-pi+aod_rx));
%             dist2 = abs((y-rx_loc_n(2))/(x-rx_loc_n(1)) - tan(rx_orient_n+aoa_rx));
%             dist3 = abs(norm([x,y]-rx_loc_n) + norm([x,y]-tx_loc) - norm(rx_loc_n-tx_loc) - (tof_rx-tof_tx)*c);

            pos_pdf(rxId, pId-1, posId, :, :) = [dist1;dist2;dist3];
        end
    end
end
% for each position and each rx, find the strongest path 
pos_feature = squeeze(min(pos_pdf,[],2)); %numRx*numPos*numPdf*G
pos_feature = permute(pos_feature,[2,4,3,1]);
pos_feature = reshape(pos_feature,size(pos_feature,1),size(pos_feature,2),[]); %numPos*numSample*nFeature
% pos_prob = prod(pos_prob,4); %numRx*1*numPos*1*G
% pos_prob = prod(pos_prob,1); %1*1*numPos*1*G
% pos_prob = squeeze(pos_prob); %numPos*G

end

