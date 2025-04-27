function psf = definePsf(N,width,bin,timeRes,c)     % refer to lct reconstruction
    linexy = linspace(-2*width,2*width,2*N-1);            % linspace of scan range
%     linexy = linspace(-2*width,2*width,2*N);    
    range = (bin*timeRes*c/2)^2;                    % total range of t^2 domain: bin^2*timeRes^2
    gridrange = bin*(timeRes*c/2)^2;                % one grid in t^2 domain: bin*timeRes^2
    [xx,yy,squarez] = meshgrid(linexy,linexy,0:gridrange:range);
    blur = abs((xx).^2+(yy).^2-squarez+0.0000001);
    blur = double(blur == repmat(min(blur,[],3),[ 1 1 bin+1 ]));
    blur = blur(:,:,1:bin);                               % generate light-cone
    psf = zeros(2*N-1,2*N-1,2*bin); 
%     psf = zeros(2*N,2*N,2*bin); 
%     psf(2:2*N,2:2*N,bin+1:2*bin) = blur;
    psf(:,:,bin+1:2*bin) = blur;                          % place it at center
end