function x = hard_threshold( x,lambda )
x(abs(x)<lambda) = 0;
end

