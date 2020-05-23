function [aoa,range] = cal_para_groundtruth(rx_loc,person_data)
% cal aoa and range from ground truth
%   rx_loc=[a;b]
%   size(person_data)=(len,2) 1 for x and 2 for y

len = size(person_data,1);
aoa=zeros(len,1);
range=zeros(len,1);
for i=1:len
    aoa(i)=(person_data(i,1)-rx_loc(1))/(person_data(i,2)-rx_loc(2));
    range(i)=sqrt((person_data(i,1)-rx_loc(1))^2+(person_data(i,2)-rx_loc(2))^2);
end

end

