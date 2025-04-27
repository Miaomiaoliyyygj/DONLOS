function [x,ii,objective] = imagedomain(d0,psf,jit,box,samp,lambdad,lambdax,mu1,mu2,mu3,au,bu,ad,bd,tau,tolerance)
    % display plots
    plots = 0;
    d0 = gpuArray(single(d0));
    jit = gpuArray(single(jit));
    
    % retrieve arguments
    [bin, N, N] = size(d0);
    max_iters = 1;
    x_iters = 1;
    d_iters = 1;
    out_iters = 1;
    mask = definemask(N,samp);
    mask = reshape(mask,[1,N,N]);
    mask = repmat(mask,[bin 1 1]);
    mask = gpuArray(single(mask)); 
    bin_resolution = 32e-12;            % Time resolution
    wall_size = 2; 
    snr = 0.1;
    width = wall_size/2;                % scan range -width to width (unit:m)
    c = 3*10^8;        % speed of light
    range = bin.*c.*bin_resolution;  %
    slope = width./range;
    lambdad = 1;
    lambdax = 1;
    z_start = 0; 
    z_stop = 1./slope;
    S = 1;
    lambda = 1;
    [mtx,mtxi] = resamplingOperator(bin);
    mtx = full(mtx);
    mtxi = full(mtxi);
    mtx = gpuArray(single(mtx));
    mtxi = gpuArray(single(mtxi));
    psf = gpuArray(single(psf)); 
    d = d0;
    vol = d0;
    
    % utility functions
    trim_array = @(x) x(S*bin/2+1:end-S*bin/2, N/2+1:end-N/2, N/2+1:end-N/2);
    pad_array = @(x) padarray(x, [S*bin/2, N/2, N/2]);
    square2cube = @(x) reshape(x, [],N,N);
    cube2square = @(x) x(:,:);
    vec = @(x) x(:);
    Ffun = @(x)  fftn(x);
    Ftfun = @(x) ifftn(x);

    % gradient kernels
    D2 = [0 0 0; 0 1 -1; 0 0 0];
    D2 = padarray(D2, [0,0,1]);
    D1 = [0 0 0; 0 1 0; 0 -1 0];
    D1 = padarray(D1, [0,0,1]);
    D3 = zeros(3,3,3);
    D3(2,2,2) = 1; D3(2,2,3) = -1;

    % operator functions
    p2o = @(x) psf2otf(x, [bin,N,N]);
    d1FT = p2o(D1);
    d2FT = p2o(D2);
    d3FT = p2o(D3);
        
    DtD = abs(d1FT).^2 + abs(d2FT).^2 + abs(d3FT).^2;
   
    pad_array1 = @(x) padarray(x, [bin/2, N/2-1, N/2-1],'pre');
    pad_array2 = @(x) padarray(x, [bin/2, N/2+1, N/2+1],'post');
    trim_array1 = @(x) x(bin/2+1:3*bin/2, N/2:3*N/2-1, N/2:3*N/2-1);
    pad_arrayx = @(x) padarray(x, [2*bin-1,2*N-1,2*N-1],'post');
    pad_arraypsf = @(x) padarray(x, [bin-1,N-1,N-1],'post');
    trim_arrayx = @(x) x(bin:2*bin-1, N:2*N-1, N:2*N-1);
    trim_arrayjit = @(x) x(4:3+bin,:,:);
    psf = permute(psf,[3 2 1]);
    psf1 = padarray(psf,[0,1,1],'pre');
    psfFT = fftn(psf1);
    psfFT1 = fftn(padarray(flip(flip(flip(psf,1),2),3),[0,1,1],'pre'));
    jit = permute(jit,[3 2 1]);
    jitFT = fftn(padarray(jit,[bin-1,N-1,N-1],'post'));
    
    A = @(x) real(trim_arrayjit(Ftfun(jitFT.*Ffun(padarray(square2cube(mtxi*cube2square(trim_array1(real(ifftshift(Ftfun(psfFT .* Ffun(pad_array2(pad_array1(x))))))))),[6 0 0],'post'))))).*mask;  
    AT = @(x) trim_array1(real(ifftshift(Ftfun(psfFT1 .* Ffun(pad_array2(pad_array1(square2cube(mtx*cube2square(x)))))))));

%     [Anorm, Asmall]= compute_operator_norm(A, AT, [bin N N]);
%     save('data/dualdomain_Anorm_HRjitbowling4.mat', 'Anorm', 'Asmall');
        
    load('data/dualdomain_Anorm_HRjitbowling4.mat');
    Anormmask = Anorm;
    fprintf('Operator eigenvalue ratio: %.02f\n', Anorm/Asmall);

    shrinkage = @(a,kappa) repmat(max(0, 1-kappa./sqrt(a(:,:,:,1).^2 + a(:,:,:,2).^2 + a(:,:,:,3).^2)),1,1,1,3).*a;

    
    % ---- For acceptance criterion ---
    alphainitmask = 2* (Anormmask).^2;
    acceptpast = 10;
    converged = 0;
    ii = 1;
    jj = 1;

    xinit = d0*0;
    x = xinit;    lambda3 = d0*0;
    lambda110 = d0*0; lambda120 = d0*0;  lambda130 = d0*0;
    lambda11 = d0*0; lambda12 = d0*0;  lambda13 = d0*0;
    lambda21 = d0*0; lambda22 = d0*0;  lambda23 = d0*0; 
    x1 = d0*0; x2 = d0*0; x3 = d0*0; d1 = d0*0; d2 = d0*0; d3 = d0*0; 
    Cu = x*0;Cd = x*0;
    % ---- For x subproblem - d initialization ---
    xprevious = x;
    Ax = A(x);
    alphamask = alphainitmask;
    Axprevious = Ax;
    grad = AT(Ax-d0);
    objective = [];
    objective = zeros(max_iters+1,1);
    objective = gpuArray(single(objective));
    objective(1) = computeobjectiveX(x,d0,d,mask,Cd,Cu,Ax,lambdax,lambdad,tau,'gaussian',1e-10,'tv',[]);
    F10 = (mu1)*DtD + alphamask;F10 = gpuArray(single(F10));
    tprevious = 1;
    yk = x;
    
    while (ii <= max_iters) && not(converged)  
        
        tic
        past = (max(ii-1-acceptpast,0):ii-1) + 1;
        maxpastobjective = max(objective(past));
        accept = 0;
        % --- Compute the step, and perform Gaussian 
        %     denoising subproblem ----
        dx = xprevious;
        z = yk - grad./alphamask;
        %compute Cu
        ux=x([2:end 1],:,:) - x(:,:,:);
        uy=x(:,[2:end 1],:) - x(:,:,:);
        uz=x(:,:,[2:end 1]) - x(:,:,:);      

        Unorm=sqrt(ux.^2+uy.^2+uz.^2);
        Unorm(Unorm==0) = 1;
        Cu = au+bu*(dxb(ux./Unorm)+dyb(uy./Unorm)+dzb(uz./Unorm)).^2; % K

        %%%%%%%%%%%%%%%%  For x

        xx = dxf(x) - lambda110/mu1;
        xy = dyf(x) - lambda120/mu1;
        xz = dzf(x) - lambda130/mu1;
        xf = sqrt(xx.^2 + xy.^2 + xz.^2);
        xf(xf==0) = 1;
        xf = max(xf - Cu/mu1,0)./xf;
        x1 = xx.*xf;
        x2 = xy.*xf;
        x3 = xz.*xf;

        g = alphamask*z - dxb(mu1*x1 + lambda110) - dyb(mu1*x2 + lambda120) - dzb(mu1*x3 + lambda130);
        g = fftn(g);
        x = real(ifftn(g./F10));
        x = max(0,x); %%%%%%%%%%%%%%

        lambda110_old = lambda110;
        lambda120_old = lambda120;
        lambda130_old = lambda130;

        lambda110 = lambda110 + mu1*(x1 - dxf(x));
        lambda120 = lambda120 + mu1*(x2 - dyf(x));
        lambda130 = lambda130 + mu1*(x3 - dzf(x));
        

        t = (1+sqrt(1+4*tprevious^2))/2;
        yk = x + ((tprevious-1)/t)*(x-xprevious);

        dx = x - dx;
        Adx = Axprevious;
        Ax = A(yk);
        Adx = Ax - Adx;
        normsqdx = sum( dx(:).^2 );

        objective(ii + 1) = computeobjectiveX(x,d0,d,mask,Cd,Cu,Ax,lambdax,lambdad,tau,'gaussian',1e-10,'tv',[]);
        
        grad = AT(Ax-d0);
        tprevious = t;
        xprevious = x;
        Axprevious = Ax; 
        converged = (abs(objective(ii + 1)-objective(ii))./abs(objective(ii+1)) <= tolerance);
        ii = ii + 1;
        
        times(ii) = toc; 
        
        vol(:,:,:) = gather(abs(x(:,:,:)));
        vol(end-50:end, :, :) = 0;
        
%         fprintf('Iteration: %d, Elapsed: %.01f\n', ii, times(ii));

        figure(1);draw3D(gather(vol),0.5,0.2,1);drawnow
        if mod(ii,30) == 0
            fprintf('Iteration: %d, Elapsed: %.01f\n', ii, times(ii));
        end
    end
    x = gather(x);
end

function objective = computeobjectiveX(x,y,d,mask,Cd,Cu,Ax,lambdax,lambdad,tau,noisetype,logepsilon,...
    penalty,varargin)
% Perhaps change to varargin 
% 1) Compute log-likelihood:
switch lower(noisetype)
    case 'poisson'
        precompute = y.*log(Ax + logepsilon);
        objective = sum(Ax(:)) - sum(precompute(:));
    case 'gaussian'
        objective = sum( (y(:) - Ax(:)).^2)./2;
    case 'gaussian2'
        objective = lambdax*sum( (d(:) - Ax(:)).^2)./2;
    case 'total'
        d1 = d.*mask;
        objective = lambdax*sum( (d(:) - Ax(:)).^2)./2+lambdad*sum( (y(:) - d1(:)).^2)./2;
end
% 2) Compute Penalty:
switch lower(penalty)
    case 'canonical'
        objective = objective + sum(abs(tau(:).*x(:)));
    case 'tv'
        x1 = dxf(x);
        x2 = dyf(x);
        x3 = dzf(x);
        objective = objective + sum(abs(Cu(:).*x1(:))) + sum(abs(Cu(:).*x2(:))) + sum(abs(Cu(:).*x3(:)));
%         objective = objective + tau.*tlv(x,'l1');
    case 'total'
        x1 = dxf(x);
        x2 = dyf(x);
        x3 = dzf(x);
        d1 = dxf(d);
        d2 = dyf(d);
        d3 = dzf(d);
        objective = objective + sum(abs(Cu(:).*x1(:))) + sum(abs(Cu(:).*x2(:))) + sum(abs(Cu(:).*x3(:)))+...;
         sum(abs(Cd(:).*d1(:))) + sum(abs(Cd(:).*d2(:))) + sum(abs(Cd(:).*d3(:)));
end
end
