clear all;
close all;
tic
%% Define parameters
N = 64*64;
width = 1;
sett = 0;% set 0 if N, M, width, bin_resolution change
snr =0.5;
z_trim = 0;%
z_offset =0;
c = 3e8;
bin_resolution =32e-12;
sampling_coeff = 1.8;

%% Load data
load('D:\mjy\mask001_data.mat');
rect_data = double(out);
%% Define PSF MTXI
N = size(rect_data,1);        % Spatial resolution of data
M = size(rect_data,3);        % Temporal resolution of data
range = M.*c.*bin_resolution; % Maximum range for histogram
savepath_psf = strcat('D:\mjy\savepsf\','psf_',num2str(N),'_',num2str(M),'_w',num2str(width),'.mat');
savepath_mtx = strcat('D:\mjy\savemtx\','mtx_',num2str(M),'.mat');
switch sett
%run if new
    case {0}
    psf = definePsf(N,M,width./range);
    [mtx,mtxi] = resamplingOperator(M);
    save(savepath_psf,'psf','-v7.3') ;
    save(savepath_mtx,'mtx', 'mtxi','-v7.3') ;
%alread exist
    case {1}
        psf_name = strcat('D:\mjy\savepsf\','psf_',num2str(N),'_',num2str(M),'_w',num2str(width),'.mat');
        mtx_name = strcat('D:\mjy\savemtx\','mtx_',num2str(M),'.mat');
        load (psf_name);
        load (mtx_name);
end

%% Reconstrcution
%% LCT
%  [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction1(rect_data,width,snr,z_trim,z_offset,psf,mtx,mtxi,bin_resolution);
%% PF
%  [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction2(rect_data,width,snr,z_trim,z_offset,psf,mtx,mtxi,bin_resolution,sampling_coeff);
%% F-K
%  [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction3(rect_data,width,z_trim,z_offset,bin_resolution);
%% FBP
% [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction4(rect_data,width,snr,z_trim,z_offset,psf,mtx,mtxi,bin_resolution);
%% SOCR
% parpool(32);
% %% adding path
% addpath('./inverse problem solver')
% addpath('./solution viewer')
% addpath('./save results')
% con_Sig= reshape(rect_data,[64*64,512]);
% con_Sig = con_Sig / max(con_Sig(:)) * 255; con_Sig = bound(con_Sig);
% %% basic settings
% basic_settings
% compute_geometry
% clc
% %% solve inverse problem
% [d,u,uDC] = full_inverse( lambda_d, lambda_pd, lambda_sd, lambda_pu,...
% loop, N_max, nped, C_i, con_Sig, kernel, pcg_iter, pcg_R_tol,...
% Arr_time,Ask_value,xmin,ymin,zmin,xg,yg,zg,Spar_init_loop,Update_u_loop,...
% k_sparse, split_Bregman_parameter, sigma_d, lambda_u_parameter,...
% lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
% ws, nno, sparse_value, sparse_threshold,...
% lambda_xyz_o, threshold_3D_o,...
% pxd, pyd, pzd, sxd, syd, szd,...
% cut_x, cut_y, cut_z);
% %% saving results
% save d d
% save u u
% %% viewing solutions
% view_solutions
% vol = permute(albedo, [3, 2, 1]);
% tic_z = linspace(0,range./2,size(vol,1));
% tic_y = linspace(-width,width,size(vol,2));
% tic_x = linspace(-width,width,size(vol,3));
%% Cur-NLOS
%  addpath('./util')
% bin = size(rect_data,3); 
% FWHM =70e-12;                      % FWHM of system temporal jitter
% sigmaBin = FWHM/bin_resolution/2/sqrt(2*log(2));
% jit = normpdf(-3:3, 0, sigmaBin);
% jit = reshape(jit, [1 1 7]);   
% y = permute(rect_data,[3 2 1]);
% grid_z = repmat(linspace(0,1,bin)',[1 N N]);
% y = y.*(grid_z.^4);
% box = drawbox(0.6,3);
% lambdad = 300; lambdax = 1; mu1 = 1; mu2 = 8e2; mu3 = 2;ad = 1e-4; bd = 1.2e-2; tau = 6e-4;bu = 3e-4;mu = 0.1;au0 = 1e-5;bu0 = 1e-5;
% tolerance = 1e-6;samp =64;
% psf = definePsf_cur(N,width,bin,bin_resolution,c);        
%  [x,iter,objectiveout] = imagedomain(y,psf,jit,box,samp,lambdad,lambdax,mu,mu2,mu3,au0,bu0,ad,bd,tau,tolerance);
% vol = x;
% tic_z = linspace(0,range./2,size(vol,1));
% tic_y = linspace(-width,width,size(vol,2));
% tic_x = linspace(-width,width,size(vol,3));
%% DO-NLOS
 [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction5(rect_data,width,z_trim,z_offset,psf,mtx,mtxi,bin_resolution);
%% View result
 vol = flip(vol, 1);
 vol = flip(vol, 3);
 set(0,'DefaultFigureVisible','on');
figure('pos',[10 10 900 300]);

subplot(1,3,1);
a1 = mat2gray(squeeze(max(vol,[],1)));
imagesc(tic_x,tic_y,a1);
title('Front view');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('x (m)');
ylabel('y (m)');
colormap('hot');
axis square;

subplot(1,3,2);
a2 = mat2gray(squeeze(max(vol,[],2)));
imagesc(tic_x,tic_z,a2);
title('Top view');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_z),max(tic_z),3));
xlabel('x (m)');
ylabel('z (m)');
colormap('hot');
axis square;

subplot(1,3,3);
a3 = mat2gray(squeeze(max(vol,[],3))');
imagesc(tic_x,tic_z,a3);
title('Side view');
set(gca,'XTick',linspace(min(tic_z),max(tic_z),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('z (m)');
ylabel('y (m)');
colormap('hot');
axis square;

toc
