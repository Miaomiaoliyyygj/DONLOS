function [Sig,con_Sig] = forward( C_i, u, N_max,...
    xg,yg,zg,xmin,ymin,zmin,Arr_time,Ask_value,kernel)
lxs = length(xg); lys = length(yg); lzs = length(zg);
u = reshape(u,lxs,lys,lzs,3);
num_measure = size(C_i,1);
Sig = NaN * zeros(num_measure,N_max);
num_grid = lxs * lys * lzs;
parfor s = 1:num_measure
       id_int = C_i(s,:);
       Arr1 = Arr_time(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1);
       Ask1 = Ask_value(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1,:);
       temp = sum( Ask1 .* u , 4);
       Sig(s,:) = full(sparse(ones(1,num_grid),Arr1(:),temp(:),1,N_max));
end
con_Sig = convolution_signal( Sig,kernel );
end