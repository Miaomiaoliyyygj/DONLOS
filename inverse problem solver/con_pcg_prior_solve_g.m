function x = con_pcg_prior_solve_g(Loss_Sig,xmin,ymin,zmin,xg,yg,zg,equi_parameter,b,x0,kernel,...
    num_measure, N_max, num_grid, C_i, Ask_value, Arr_time, max_iter, max_r_tol)
lxs = length(xg); lys = length(yg); lzs = length(zg); l_ker = length(kernel);

nb = norm(b); r_err = NaN * zeros(1,max_iter);
x = x0;
r = b - AA( x,equi_parameter,...
num_measure,num_grid,kernel,l_ker,Ask_value,...
zmin,zg,ymin,yg,xmin,xg,Arr_time,C_i,lzs,lys,lxs,N_max );
rr = r' * r; p = r; r_err(1) = norm(r)/nb;
if r_err(1) >= max_r_tol
    for k = 2:max_iter   
        AAp = AA( p,equi_parameter,...
num_measure,num_grid,kernel,l_ker,Ask_value,...
zmin,zg,ymin,yg,xmin,xg,Arr_time,C_i,lzs,lys,lxs,N_max );
        pAAp = p' * AAp;
        alpha = rr / pAAp;
        x_new = x + alpha * p;
        r_new = r - alpha * AAp;        
        r_err(k) =  norm(r_new)/nb;
        if r_err(k) < max_r_tol
            x = x_new;
            break
        end        
        rr_new = r_new' * r_new;
        beita = rr_new / rr;
        p = r_new + beita * p;
        x = x_new;
        r = r_new;
        rr = rr_new;
    end
end

function y = AA( x,equi_parameter,num_measure,num_grid,kernel,l_ker,Ask_value,zmin,zg,ymin,yg,xmin,xg,Arr_time,C_i,lzs,lys,lxs,N_max )
x = reshape(x,lxs,lys,lzs,3);
y = zeros(lxs,lys,lzs,3);    
parfor s = 1:num_measure
    if Loss_Sig(s) == 1, continue, end
    mid = zeros(N_max,1);
    temp_y = zeros(lxs,lys,lzs,3); id_int = C_i(s,:);
    Arr1 = Arr_time(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1);
    Ask1 = Ask_value(xg - id_int(1) - xmin + 1, yg - id_int(2) - ymin + 1, zg - id_int(3) - zmin + 1,:);        
    for conv_bin = 1:l_ker        
            temp = sum( kernel(conv_bin) * Ask1 .* x, 4); temp = temp(:);
            mid = full(sparse( Arr1(:) + conv_bin - 1,ones(num_grid,1),temp, N_max, 1));
    end
    for conv_bin = 1:l_ker
    temp_y = kernel(conv_bin) * Ask1(:,:,:,:) .* repmat(mid(Arr1 + conv_bin - 1),1,1,1,3);
    end
    y = y + temp_y;
end 
y = equi_parameter * x + y;
y = y(:);
end

end