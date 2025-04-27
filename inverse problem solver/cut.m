function u = cut( u, s1, s2, s3 )
u(1:s1,:,:) = 0; u(:,1:s2,:) = 0; u(:,:,1:s3) = 0;
u(end-s1+1:end,:,:) = 0; u(:,end-s2+1:end,:) = 0; u(:,:,end-s3+1:end) = 0;
end

