function [d,u,uDC] = full_inverse( lambda_d, lambda_pd, lambda_sd, lambda_pu,...
loop, N_max, nped, C_i, con_Sig, kernel, pcg_iter, pcg_R_tol,...
Arr_time,Ask_value,xmin,ymin,zmin,xg,yg,zg,Spar_init_loop,Update_u_loop,...
k_sparse, split_Bregman_parameter, sigma_d, lambda_u_parameter,...
lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
ws, nno, sparse_value, sparse_threshold,...
lambda_xyz_o, threshold_3D_o,...
pxd, pyd, pzd, sxd, syd, szd,...
cut_x, cut_y, cut_z)
%% Step 0: Signal
if size(con_Sig,2) < N_max
    con_Sig = [con_Sig,zeros(size(con_Sig,1),N_max - size(con_Sig,2))];
end
total_time = 0;
Loss_Sig = ( sum(con_Sig,2) == 0 );
%% Step 1: Parameters
lxs = length(xg); lys = length(yg); lzs = length(zg);
[num_measure,N_sig] = size(con_Sig); num_grid = lxs * lys * lzs;
%% Step 2: Initialization d
d = con_Sig;
%% Step 3: Initialization u
fprintf('----------------------------------------\n')
fprintf('u = argmin_u ||Au-b||^2\n')
tic
RHS = generate_RHS_g(Loss_Sig,xmin,ymin,zmin,xg,yg,zg,C_i,Arr_time,Ask_value,d,kernel);
u = con_pcg_solve_g(Loss_Sig,xmin,ymin,zmin,xg,yg,zg,kernel,num_measure, N_max,num_grid, C_i, Ask_value, Arr_time,zeros(lxs*lys*lzs*3,1), pcg_iter, RHS,pcg_R_tol);
toc; 
total_time = total_time + toc;
u = reshape(u,num_grid,3);
cd('save results'), save('u00.mat','u'), cd ..
%% Step 4: Forward operator
fprintf('----------------------------------------\n')
fprintf('computing Au\n')
tic
[~,d_u] = forward( C_i, u, N_sig, xg,yg,zg,xmin,ymin,zmin,Arr_time,Ask_value,kernel);
toc
total_time = total_time + toc;
cd('save results'), save('d0.mat','d_u'), cd ..
%% Step 5: Choosing parameters
fprintf('----------------------------------------\n')
fprintf('choosing parameters\n')
data_error = sum(sum((d - d_u).^2)); data_L21_norm = sum(row_norm( u ));
s_u = k_sparse * data_error / data_L21_norm;
% choose mu
row_norm_u = row_norm( u ); row_norm_u(row_norm_u == 0) = [];
mu = split_Bregman_parameter * s_u * mean(1 ./ row_norm_u );
clear row_norm_u
fprintf('Data Misfit : %1.3g\n', sqrt(data_error))
fprintf('u_{2,1} norm : %1.3g\n', data_L21_norm)
fprintf('s_u chosen : %1.3g\n', s_u)
fprintf('mu chosen : %1.3g\n', mu)
%% Step 6: Split Bregman Iteration for a sparse solution
fprintf('----------------------------------------\n')
fprintf('u = argmin_u ||Au-b||^2 + ||u||_{2,1}\n')
b = zeros(num_grid,3);
R_uv = NaN * zeros(1,Spar_init_loop);
tic
for l1_loop = 1:Spar_init_loop
% update v
temp_var = u - b;
v = repmat(max(0,1 - s_u / mu / 2 ./row_norm(temp_var)),1,3) .* temp_var;
% update u
temp_var = b + v; temp_var = temp_var(:);
u = con_pcg_prior_solve_g(Loss_Sig, xmin,ymin,zmin,xg,yg,zg,mu,RHS + mu * temp_var,v(:),kernel,num_measure, N_max, num_grid, C_i, Ask_value, Arr_time, pcg_iter, pcg_R_tol);
u = reshape(u,num_grid,3);
R_uv(l1_loop) = norm(u(:) - v(:)) / norm(u(:));
cd('save results'), save('Ruv0.mat','R_uv'), cd ..
if R_uv(l1_loop) < 5E-3
fprintf('At iteration %d, R(u,v) = %1.2f\n', l1_loop, R_uv(l1_loop));
    break
end
% update b
b = b + v - u;
end
toc
total_time = total_time + toc;
u = v;
cd('save results'), save('u0.mat','u'), cd ..
%% Step 7: Initialization dictionaries
u = reshape(u,lxs,lys,lzs,3);
temp = u(:,:,:,1); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,1) = temp;
temp = u(:,:,:,2); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,2) = temp;
temp = u(:,:,:,3); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,3) = temp; clear temp
u = reshape(u,[],3);
tic
denoise_factor = 255 / max( row_norm(u) );
fprintf('Initializing dictionaries\n')
[Dz,Dxy,uDC] = denoise_albedo( u, lxs,lys,lzs, denoise_factor, lambda_pu,...
                               lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
                               ws, nno, sparse_value, sparse_threshold,...
                               lambda_xyz_o, threshold_3D_o);
toc
total_time = total_time + toc;
cd('save results'), save('Dz0.mat','Dz'), cd ..
cd('save results'), save('Dxy0.mat','Dxy'), cd ..
cd('save results'), save('uDC0.mat','uDC'), cd ..
%% Step 8: Choosing parameters
denoise_u_error = sum(sum( (u - uDC).^2 ));
lambda_u = lambda_u_parameter * data_error / denoise_u_error;
fprintf('||Au - d|| : %1.3g\n', sqrt(data_error))
fprintf('||u - uDC|| : %1.3g\n', sqrt(denoise_u_error))
fprintf('lambda_u : %1.3g\n', lambda_u)
%% Step 9: Main loop
for num_iter = 1:loop
      num_iter
% update d
tic
fprintf('update d\n')
[~,d_u] = forward( C_i, u, N_sig,xg,yg,zg,xmin,ymin,zmin,Arr_time,Ask_value,kernel);
d = update_d( con_Sig, d_u, lambda_d, lambda_pd, lambda_sd, nped, N_sig, N_max, sigma_d, pxd, pyd, pzd, sxd, syd, szd );
toc
total_time = total_time + toc;
cd('save results'), save(['d',num2str(num_iter),'.mat'],'d'), cd ..
% update u
tic
u = reshape(u,num_grid,3);
fprintf('update u\n');
b = zeros(num_grid,3);
RHS = generate_RHS_g( Loss_Sig,xmin,ymin,zmin,xg,yg,zg,C_i,Arr_time, Ask_value,d,kernel );
R_uv = NaN * zeros(1,Update_u_loop);
for upu_loop = 1:Update_u_loop 
% update v
temp_var = u - b;
v = repmat(max(0,1 - s_u / mu / 2 ./row_norm(temp_var)),1,3) .* temp_var;
% update u
temp_var = b + v;
temp_var = temp_var(:);
u = con_pcg_prior_solve_g(Loss_Sig,xmin,ymin,zmin,xg,yg,zg,lambda_u + mu,RHS + lambda_u * uDC(:) + mu * temp_var, ...
    v(:),kernel, num_measure, N_max, num_grid, C_i, Ask_value, Arr_time, pcg_iter, pcg_R_tol);
u = reshape(u,num_grid,3);
R_uv(upu_loop) = norm(u-v)/norm(u);
if R_uv(upu_loop) < 5E-3
fprintf('At iteration %d, R(u,v) = %1.2f\n', upu_loop, R_uv(upu_loop));
    break
end
% update b
b = b + v - u;
end
toc
total_time = total_time + toc;
u = v;
cd('save results'), save(['u',num2str(num_iter),'.mat'],'u'), cd ..
% update dictionaries
u = reshape(u,lxs,lys,lzs,3);
temp = u(:,:,:,1); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,1) = temp;
temp = u(:,:,:,2); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,2) = temp;
temp = u(:,:,:,3); temp = cut(temp,cut_x,cut_y,cut_z); u(:,:,:,3) = temp; clear temp
u = reshape(u,[],3);
tic
fprintf('update dictionaries and coefficients\n')
[Dz,Dxy,uDC] = denoise_albedo( u, lxs,lys,lzs, denoise_factor, lambda_pu,...
                               lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
                               ws, nno, sparse_value, sparse_threshold,...
                               lambda_xyz_o, threshold_3D_o);
toc
total_time = total_time + toc;
cd('save results'), save(['Dz',num2str(num_iter),'.mat'],'Dz'), cd ..
cd('save results'), save(['Dxy',num2str(num_iter),'.mat'],'Dxy'), cd ..
cd('save results'), save(['uDC',num2str(num_iter),'.mat'],'uDC'), cd ..
end
fprintf('----------------------------------------\n')
fprintf('total time = %1.3g s\n',total_time)
end