function noisy_img = Wiener_sig_quick( basic_img, noisy_img,lambda_sd, sigma_d,...
                                       pxd, pyd, pzd, sxd, syd, szd)
sigma_d = sigma_d / (sqrt(1 + lambda_sd));
[m,n,o] = size(noisy_img);
mpx = m - pxd + 1; npy = n - pyd + 1; opz = o - pzd + 1;
Px = 1:sxd:mpx; if Px(length(Px))~=mpx, Px = [Px,mpx]; end, lx = length(Px);
Py = 1:syd:npy; if Py(length(Py))~=npy, Py = [Py,npy]; end, ly = length(Py);
Pz = 1:szd:opz; if Pz(length(Pz))~=opz, Pz = [Pz,opz]; end, lz = length(Pz);
Dxy = kron(dctmtx(pzd) , kron(dctmtx(pyd) , dctmtx(pxd)))';
Pun = NaN * zeros(pxd * pyd * pzd, lx * ly * lz);
Pub = NaN * zeros(pxd * pyd * pzd, lx * ly * lz);
row = 0;
count_num = zeros(m,n,o);
for j = 0:pxd-1
    for k = 0:pyd-1
        for l = 0:pzd-1
            row = row + 1;
            temp_noisy = noisy_img(Px+j,Py+k,Pz+l);
            temp_basic = basic_img(Px+j,Py+k,Pz+l);
            count_num(Px+j,Py+k,Pz+l) = count_num(Px+j,Py+k,Pz+l) + 1;
            Pun(row,:) = temp_noisy(:)';
            Pub(row,:) = temp_basic(:)';
        end
    end
end
basic_coef = Dxy' * Pub; 
noisy_coef = Dxy' * Pun;
noisy_coef = ( basic_coef.^2 ./ (basic_coef.^2 + sigma_d^2) ) .* noisy_coef;
Pun = Dxy * noisy_coef;
row = 0;
for j = 0:pxd-1
    for k = 0:pyd-1
        for l = 0:pzd-1
            row = row + 1;
            noisy_img(Px+j,Py+k,Pz+l) = noisy_img(Px+j,Py+k,Pz+l) + reshape(Pun(row,:),lx,ly,lz); 
        end
    end
end
noisy_img = noisy_img ./ count_num;
noisy_img(noisy_img < 0) = 0;
end