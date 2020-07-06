%% creates set of steering vectors. MUSIC spectrum is calculated at each of these steering vectors.
function [aTot] = xyVecBackscatter3D(antenna_type,calib_const,generateAtot)
% % % Input
% cc.x_range: range of x axis
% cc.x_step: gap in x axis
% cc.y_range: range of y axis
% cc.y_step: gap in y axis
% cc.AS_vec: position of rx antennas
% cc.DS_vec: position of rx antennas
% generateAtot: whether we need to generate a new Atot or load the previous one
% calib_const: parameters
% % % Output
% aTot: matrix with each column being a steering vector (nSubcarr*nRx,nX*nY);
cc = calib_const;
nX = floor((cc.x_range(2)-cc.x_range(1))/cc.x_step);
nY = floor((cc.y_range(2)-cc.y_range(1))/cc.y_step);

rAPs = {'rAP1','rAP2','rAP3'};
if any(strcmp(rAPs,antenna_type)) % calculate AoA
    APid = str2num(antenna_type(4));
    nTx = size(cc.DS_vec,2); % number of antennas on Tx
    nRx = size(cc.AS_vec(APid,:),2); % number of antennas on Rx

    nSubCarr = cc.F;

    aTot = zeros(nSubCarr*nRx,nX*nY);
    for yi = 1:nY
        for xi=1:nX
            x = cc.x_range(1)+(xi-1/2)*cc.x_step;
            y = cc.y_range(1)+(yi-1/2)*cc.y_step;
            aoa = atan2(y-cc.rx_loc(APid,2), x-cc.rx_loc(APid,1))-(cc.rx_orient(APid)-pi);
            delay = norm([x,y]-cc.tx_loc)+norm([x,y]-cc.rx_loc(APid,:));
            steeringVec = gridSampleBackscatter3D(delay,aoa,0,cc.AS_vec(APid,:),cc.DS_vec,calib_const);
            id = (yi-1)*nX+xi;
            aTot(:,id)=steeringVec;
        end
    end
end

end

