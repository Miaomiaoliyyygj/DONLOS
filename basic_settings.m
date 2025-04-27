%% basic parameters
ext_x = 0; ext_y = 7; ext_z = 7;
cut_x = 0; cut_y = 0; cut_z = 0;
show_cut_x = 7; show_cut_y = 7; show_cut_z = 7;
c = 3E8; delta_t = 32E-12;
bx = c * delta_t; by = 2/63; bz = 2/63;
xb = 1.20; xe = 1.50; yb = 0; ye = 2; zb = 0; ze = 2;
nbx = 1; nby = 1; nbz = 1;
nped = 64;
niy = 1; ndy = 1;
niz = 1; ndz = 1;

%% reconstruction parameters
k_sparse = 2.5;
split_Bregman_parameter = 0.5;
lambda_u_parameter = 1;
lambda_pu = 40;
sigma_d = 80;
lambda_d = 1;
lambda_sd = 1/4;
lambda_pd = 16;
loop = 1;
Spar_init_loop = 20;
Update_u_loop = 3;
pcg_iter = 20;
pcg_R_tol = 5E-3;
%% additional parameters
lambda_RMSE = 50;                                                          
pxo = 1; pyo = 2; pzo = 2;                                                    
sxo = 1; syo = 1; szo = 1;                                                    
ws = 6;                                                  
nno = 3;                                                                  
sparse_value = 40;                                                         
sparse_threshold = 2; 
lambda_xyz_o = 2.6; 
threshold_3D_o = 2.6; 
pxd = 3; pyd = 3; pzd = 3;                                                    
sxd = 1; syd = 1; szd = 1;
%% illumination points
ixbc = 0; ixec = 0;
iybc = 0; iyec = (nped - 1) * niy;
izbc = 0; izec = (nped - 1) * niz;
C_i = NaN * zeros(nped^2,3);
col = 0;
for y_axis = 1:nped
    for z_axis = 1:nped
        col = col + 1;
        C_i(col,:) = [0, (nped - y_axis) * niz, (nped - z_axis) * niy];
    end
end
clear y_axis z_axis col
%% detection points
C_d = C_i;
%% reconstruction domain
xb = xb - bx * ext_x; yb = yb - by * ext_y; zb = zb - bz * ext_z;
xe = xe + bx * ext_x; ye = ye + by * ext_y; ze = ze + bz * ext_z;
xbc = floor( xb / bx ); xec = ceil(xe / bx); xg = xbc:nbx:xec; lxs = length(xg);
ybc = floor( yb / by ); yec = ceil(ye / by); yg = ybc:nby:yec; lys = length(yg);
zbc = floor( zb / bz ); zec = ceil(ze / bz); zg = zbc:nbz:zec; lzs = length(zg);
num_grid = lxs * lys * lzs;
%% convolution kernel
kernel = 1;
%% normal of the visible wall
nPd = [1 0 0]';
%% number of grids in each direction
xmin = xbc - ixec; xmax = xec - ixbc; numx = xmax - xmin + 1;
ymin = ybc - iyec; ymax = yec - iybc; numy = ymax - ymin + 1;
zmin = zbc - izec; zmax = zec - izbc; numz = zmax - zmin + 1;
num_geo = numx * numy * numz;