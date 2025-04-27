function Sig2D = sig3dto2d( Sig3D, nped,bin2D )
Sig2D = zeros(nped.^2,bin2D);
[~,~,bin3D] = size(Sig3D);
col = 0;
for j = 1:nped
    for k = 1:nped
        col = col + 1;
        temp_Sig = Sig3D(k,j,:);
        temp_Sig = temp_Sig(:);
        Sig2D(col,1:bin3D) = temp_Sig;      
    end
end
end