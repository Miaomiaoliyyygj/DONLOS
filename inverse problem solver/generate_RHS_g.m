function RHS = generate_RHS_g( Loss_Sig,xmin,ymin,zmin,xg,yg,zg,C_i,Arr_time, Ask_value,con_Sig,kernel )
lxs = length(xg); lys = length(yg); lzs = length(zg);
num_measure = size(con_Sig,1);
l_ker = length(kernel);
RHS = zeros(lxs,lys,lzs,3);
parfor s = 1:num_measure
       if Loss_Sig(s) == 1, continue, end
       temp_sig = con_Sig(s,:); id_int = C_i(s,:);                    
       Arr1 = Arr_time(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1);
       Ask1 = Ask_value(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1,:);
       temp_RHS = zeros(lxs,lys,lzs,3,l_ker);
       for conv_bin = 1 : l_ker 
       temp_RHS(:,:,:,:,conv_bin) = kernel(conv_bin) *  repmat(temp_sig(Arr1 + conv_bin-1),1,1,1,3);
       end
       temp_RHS = sum(temp_RHS,5) .* Ask1;    
       RHS = RHS + temp_RHS;
end
RHS = RHS(:);
end