function u = bound( u )
u(u>255) = 255;
u(u<0) = 0;
end

