
function [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction3(rect_data,width,z_trim,z_offset,bin_resolution)



% Constants
c              = 3e8;   % Speed of light (meters per second)

% Adjustable parameters

K          = 0;         % Downsample data to (4 ps) * 2^K = 16 ps for K = 2


% Load scene & set visualization parameter

N = size(rect_data,1);        % Spatial resolution of data
M = size(rect_data,3);        % Temporal resolution of data
range = M.*c.*bin_resolution; % Maximum range for histogram
    
% Downsample data to 16 picoseconds
for k = 1:K
    M = M./2;
    bin_resolution = 2*bin_resolution;
    rect_data = rect_data(:,:,1:2:end) + rect_data(:,:,2:2:end);
    z_trim = round(z_trim./2);
    z_offset = round(z_offset./2);
end
    
% Set first group of histogram bins to zero (to remove direct component)
 rect_data(:,:,1:z_trim) = 0;
 

% Permute data dimensions
data = permute(rect_data,[3 2 1]);

% Define volume representing voxel distance from wall
grid_z = repmat(linspace(0,1,M)',[1 N N]);
% for i=1:128
%     f(:,:,i)=grid_z(i,:,:);
% end
% g= zeros(256,256)
% g=f(:,:,50);
display('Inverting...');
% tic;


 % Step 0: define virtual wavelet properties
         [z,y,x] = ndgrid(-M:M-1,-N:N-1,-N:N-1);
        z = z./M; y = y./N; x = x./N;

        % Step 0: Pad data
        data = data .* grid_z.^2;
        data = sqrt(data);
        tdata = zeros(2.*M,2.*N,2.*N);
        tdata(1:end./2,1:end./2,1:end./2) = data;

        % Step 1: FFT
        tdata = fftshift(fftn(tdata));

        % Step 2: Stolt trick
        tvol = interpn(z,y,x,tdata,sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),y,x,'linear',0);
        tvol = tvol.*(z > 0);
        tvol = tvol.*abs(z)./max(sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(x.^2+y.^2)+z.^2)),1e-6);

        % Step 3: IFFT
        tvol = ifftn(ifftshift(tvol));
        tvol = abs(tvol).^2;    
        vol = abs(tvol(1:end./2,1:end./2,1:end./2));



display('... done.');
time_elapsed = toc;

display(sprintf(['Reconstructed volume of size %d x %d x %d '...
    'in %f seconds'], size(vol,3),size(vol,2),size(vol,1),time_elapsed));

tic_z = linspace(0,range./2,size(vol,1));
tic_y = linspace(-width,width,size(vol,2));
tic_x = linspace(-width,width,size(vol,3));

% Crop and flip reconstructed volume for visualization
ind = round(M.*2.*width./(range./2));
vol = vol(:,:,end:-1:1);
vol = vol((1:ind)+z_offset,:,:);
% 

tic_z = tic_z((1:ind)+z_offset);




