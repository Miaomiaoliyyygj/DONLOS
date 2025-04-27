[JJ,KK,LL] = meshgrid(ymin:ymax,xmin:xmax,zmin:zmax);
JJ = by * JJ; KK = bx * KK; LL = bz * LL;
n_vec = sqrt(JJ.^2 + KK.^2 + LL.^2);
Arr_time = ceil(2 * n_vec /c/delta_t);
Ask_value = repmat( 1 ./ (n_vec.^5) ,1,1,1,3);
Ask_value(:,:,:,1) = KK .* Ask_value(:,:,:,1);
Ask_value(:,:,:,2) = JJ .* Ask_value(:,:,:,2);
Ask_value(:,:,:,3) = LL .* Ask_value(:,:,:,3);
Ask_value = - Ask_value;
N_max = max(Arr_time(:));
clear JJ KK LL n_vec