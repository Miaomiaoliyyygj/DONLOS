function d = update_d( con_Sig, d_u, lambda_d, lambda_pd, lambda_sd, nped, bin3D, N_max, sigma_d, pxd, pyd, pzd, sxd, syd, szd )
d = (d_u + lambda_d * con_Sig) / (1 + lambda_d);
for Cd_loop = 1:2
noisy_img = (con_Sig + lambda_sd * d) / (1 + lambda_sd);
noisy_img = sig2dto3d( noisy_img, nped,bin3D );
basic_img = sig2dto3d( d_u, nped,bin3D );
denoised_img = Wiener_sig_quick( basic_img, noisy_img,lambda_sd, sigma_d, pxd, pyd, pzd, sxd, syd, szd );
denoised_img = sig3dto2d( denoised_img, nped, N_max );
d = (d_u + lambda_d * con_Sig + lambda_d * lambda_pd * lambda_sd * denoised_img) /... 
(1 + lambda_d + lambda_d * lambda_pd * lambda_sd );
end
end