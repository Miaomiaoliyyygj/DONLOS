function Y = sub_transpose( Bu,nn2D )
[pxy,Nnn2D] = size(Bu);                                                 
N = Nnn2D / nn2D;                                                         
Y = NaN * zeros(nn2D,pxy * N);                                             
for j = 1:N                                                                
    temp_block = Bu(:, nn2D*(j-1) + 1 : nn2D*(j-1) + nn2D);               
    Y(:,pxy*(j-1) + 1 : pxy*(j-1) + pxy) = temp_block';                    
end
end