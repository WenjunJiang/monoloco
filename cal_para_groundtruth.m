function [aoa,aod,range] = cal_para_groundtruth(tx_loc, rx_loc,person_data)
% cal aoa and range from ground truth
%   rx_loc=[a;b]
%   size(person_data)=(len,2) 1 for x and 2 for y

len = size(person_data,1);
numRx = size(rx_loc,1);
aoa=zeros(len,numRx);
range=zeros(len,numRx);
aod = zeros(len,1);
for i=1:len
    aod(i)=atan2(person_data(i,2)-tx_loc(2), person_data(i,1)-tx_loc(1));
    for rx=1:numRx
        aoa(i,rx)=atan2(person_data(i,2)-rx_loc(rx,2), person_data(i,1)-rx_loc(rx,1));
        range(i,rx)=sqrt((person_data(i,1)-rx_loc(rx,1))^2+(person_data(i,2)-rx_loc(rx,2))^2);
    end
end

end

