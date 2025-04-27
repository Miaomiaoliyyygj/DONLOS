function u=ITV_ROF(f,mu)
% Isotropic ROF denoising
% This function performs the minimization of
% u=arg min sqrt(|D_x u|^2+|D_y u|^2)+0.5*mu*||u-f||_2^2
% by the Split Bregman Iteration
% f = noisy image
% mu = regularization parameter
% Reference£ºGoldstein T , Osher S . The Split Bregman Method for
% L1-Regularized Problems[J].
%SIAM Journal on Imaging Sciences, 2009, 2(2):323-343.
%============================================
[M,N]=size(f);
f=double(f);
dx=zeros(M,N);
dy=zeros(M,N);
bx=zeros(M,N);
by=zeros(M,N);
s=zeros(M,N);
u=f;
Z=zeros(M,N);
lambda=100;
Mask=zeros(M,N);
Mask(1,1) = 1;
FMask=fft2(Mask);
%Fourier Laplacian mask initialization
D = zeros(M,N);
D([end,1,2],[end,1,2]) = [0,1,0;1,-4,1;0,1,0];
FD=fft2(D);
%Fourier constant initialization
FW=((mu/lambda)*abs(FMask).^2-real(FD)).^-1;
FF=(mu/lambda)*conj(FMask).*fft2(f);
err=1;
tol=1e-3;
while err>tol
    %while K<Niter,
    tx=dx-bx;
    ty=dy-by;
    up=u;
    %Update u
    u=real(ifft2(FW.*(FF-fft2(tx-tx(:,[1,1:N-1])+ty-ty([1,1:M-1],:)))));
    ux=u-u(:,[1,1:N-1]);
    uy=u-u([1,1:M-1],:);
    tmpx=ux+bx;
    tmpy=uy+by;
    s=sqrt(tmpx.^2+tmpy.^2);
    thresh=max(Z,s-1/lambda)./max(1e-12,s);
    %Update dx
    dx=thresh.*tmpx;
    %Update dy
    dy=thresh.*tmpy;
    %Update bx and by
    bx=tmpx-dx;
    by=tmpy-dy;
    err=sum(sum((up-u).^2));
end
