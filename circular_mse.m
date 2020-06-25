function dist = circular_mse(mat1,mat2)
    tmp = mod(mat1+3*pi-mat2,2*pi)-pi; % make the first item larger than the second
    dist = mean(tmp.^2,'all');
end

