function [Dz,Dxy,denoised_img] = eat_vec( noisy_img, denoise_factor, lambda_pu,...
                                          lambda_RMSE, pxo, pyo, pzo, sxo,syo, szo,...
                                          ws, nno, sparse_value, sparse_threshold,...
                                          lambda_xyz_o, threshold_3D_o)
wl = 2 * ws + 1;
lambda_xyz_o = lambda_xyz_o * lambda_pu;
threshold_3D_o = threshold_3D_o * lambda_pu; 
[m,n,o] = size(noisy_img);
noisy_img = noisy_img * denoise_factor;
N_iter = 10; N_xy_iter = 5; N_z_iter = 2;                                                     
pxyz = pxo * pyo * pzo;
mpx = m - pxo + 1; npy = n - pyo + 1; opz = o - pzo + 1;
Px = 1:sxo:mpx; if Px(length(Px))~=mpx, Px = [Px,mpx]; end, lx = length(Px);
Py = 1:syo:npy; if Py(length(Py))~=npy, Py = [Py,npy]; end, ly = length(Py);
Pz = 1:szo:opz; if Pz(length(Pz))~=opz, Pz = [Pz,opz]; end, lz = length(Pz);
Bu = NaN * zeros(pxyz,lx*ly*lz*nno);
parti_block = 0;
for j = 1:lx
    for k = 1:ly
        for l = 1:lz
            jj = Px(j); kk = Py(k); ll = Pz(l);
            ref_patch = noisy_img(jj:jj+pxo-1,kk:kk+pyo-1,ll:ll+pzo-1);
            ref_patch = ref_patch(:);
            if sum(abs(ref_patch) > sparse_value) < sparse_threshold
                continue
            end            
            count_num = 0;
            temp_block = NaN * zeros(pxyz,wl^3);
            MSE = NaN * zeros(1,wl^3);
            for x = -ws:ws
                for y = -ws:ws
                    for z = -ws:ws
                        tj = jj + x; tk = kk + y; tl = ll + z;
                        if tj >= 1 && tj <= mpx &&...        
                           tk >= 1 && tk <= npy &&...
                           tl >= 1 && tl <= opz
                           count_num = count_num + 1;
nei_patch = noisy_img(tj:tj+pxo-1,tk:tk+pyo-1,tl:tl+pzo-1);
if sum(abs(ref_patch) > sparse_value) < sparse_threshold
    continue
end         
nei_patch = nei_patch(:);
temp_block(:,count_num) = nei_patch;
MSE(count_num) = sum((ref_patch - nei_patch).^2)/pxyz;
                        end
                    end
                end
            end
 [MSE_ascend_value,MSE_ascend_order] = sort(MSE,'ascend');        
 f = find(MSE_ascend_value < lambda_RMSE.^2);
 temp_block = temp_block(:,MSE_ascend_order(1:nno));
 if length(f) < nno 
     continue
 else
    Bu(:,parti_block + 1 : parti_block + nno) = temp_block;
    parti_block = parti_block + nno;
 end
        end
    end
end
Bu(:,isnan(Bu(1,:)) == 1) = [];
Dxy = kron(dctmtx(pyo)',dctmtx(pxo)');
Dxy = kron(dctmtx(pzo)',Dxy);
Dz = dctmtx(nno)';                 
for j = 1:N_iter
    Coef_xy = sub_transpose(Dxy' * Bu,nno);
    for j_sub = 1:N_z_iter
    C = hard_threshold(Dz' * Coef_xy, lambda_xyz_o);
    [U,~,V] = svd(Coef_xy * C'); Dz = U * V';
    end
    Coef_z = sub_transpose(Dz' * sub_transpose(Bu,nno),pxyz);
    for j_sub = 1:N_xy_iter
    C = hard_threshold(Dxy' * Coef_z, lambda_xyz_o);
    [U,~,V] = svd(Coef_z * C'); Dxy = U * V';
    end      
end
ebuff = zeros(m,n,o);
wbuff = zeros(m,n,o);
for j = 1:lx
    parfor k = 1:ly
        for l = 1:lz
            jj = Px(j); kk = Py(k); ll = Pz(l);
            ref_patch = noisy_img(jj:jj+pxo-1,kk:kk+pyo-1,ll:ll+pzo-1);
            ref_patch = ref_patch(:);     
            if sum(abs(ref_patch) > sparse_value) < sparse_threshold
                continue
            end
            temp_est_mtx = zeros(m,n,o);
            temp_w_mtx = zeros(m,n,o);
            count_num = 0;
            temp_block = NaN * zeros(pxyz,wl^3);
            MSE = NaN * zeros(1,wl^3);
            nei_x = NaN * zeros(1,wl^3);
            nei_y = NaN * zeros(1,wl^3);
            nei_z = NaN * zeros(1,wl^3);
            for x = -ws:ws
                for y = -ws:ws
                    for z = -ws:ws
                        tj = jj + x; tk = kk + y; tl = ll + z;
                        if tj >= 1 && tj <= mpx &&...        
                           tk >= 1 && tk <= npy &&...
                           tl >= 1 && tl <= opz
                           count_num = count_num + 1;
nei_patch = noisy_img(tj:tj+pxo-1,tk:tk+pyo-1,tl:tl+pzo-1);

if sum(abs(ref_patch) > sparse_value) < sparse_threshold
    continue
end

nei_patch = nei_patch(:);
temp_block(:,count_num) = nei_patch;
MSE(count_num) = sum((ref_patch - nei_patch).^2)/pxyz;
nei_x(count_num) = x; nei_y(count_num) = y; nei_z(count_num) = z;
                        end
                    end
                end
            end
 [MSE_ascend_value,MSE_ascend_order] = sort(MSE,'ascend');       
 f = find(MSE_ascend_value < lambda_RMSE.^2);
 if length(f) < nno
    temp_Dz = dctmtx(length(f))';
 else
    temp_Dz = Dz;
    f = f(1:nno);
 end
 count_num = length(f);
 noisy_block = temp_block(:,MSE_ascend_order(f));
 nei_x = nei_x(MSE_ascend_order(f));
 nei_y = nei_y(MSE_ascend_order(f));
 nei_z = nei_z(MSE_ascend_order(f)); 
 temp_coef = Dxy' * noisy_block;
 temp_coef = temp_Dz' * temp_coef';
 temp_coef(abs(temp_coef) < threshold_3D_o) = 0;
 temp_weight = max(1,sum(sum(temp_coef~=0))) * count_num ;
 temp_coef = temp_Dz * temp_coef;
 denoised_block = Dxy * temp_coef';
for nei = 1:count_num
            temp_patch = reshape(denoised_block(:,nei),pxo,pyo,pzo);
            tnx = jj + nei_x(nei); tny = kk + nei_y(nei); tnz = ll + nei_z(nei);
            temp_w_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) =...
                temp_w_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) + temp_weight;
            temp_est_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) =...
                temp_est_mtx(tnx:tnx+pxo-1,tny:tny+pyo-1,tnz:tnz+pzo-1) +  temp_weight * temp_patch;
end
        ebuff = ebuff + temp_est_mtx;
        wbuff = wbuff + temp_w_mtx; 
        end
    end
end
denoised_img = ebuff./wbuff;
denoised_img(isnan(denoised_img)==1) = 0;
denoised_img(denoised_img < 0) = 0;
denoised_img = denoised_img / denoise_factor;
end