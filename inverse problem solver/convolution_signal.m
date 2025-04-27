function y = convolution_signal( Sig,kernel )
[num_measure,N_max] = size(Sig);
y = zeros(num_measure,N_max);
l_ker = length(kernel);
for j = 1:l_ker
    y = y + kernel(j) * Sig;
    Sig = Sig(:,[N_max,1:N_max - 1]);
    Sig(:,1) = 0;
end
end