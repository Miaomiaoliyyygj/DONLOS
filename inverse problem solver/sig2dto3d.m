function Sig3D = sig2dto3d( Sig2D, nped,bin3D )
Sig3D = NaN * zeros(nped,nped,bin3D);
col = 0;
for j = 1:nped
    for k = 1:nped
        col = col + 1;
        temp_Sig = Sig2D(col,:);
        temp_Sig = temp_Sig(1:bin3D);
        Sig3D(k,j,:) = temp_Sig;
    end
end
end