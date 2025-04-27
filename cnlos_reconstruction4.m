function [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction4(rect_data,width,snr,z_trim,z_offset,psf,mtx,mtxi,bin_resolution)



% Constants
c              = 3e8;   % Speed of light (meters per second)

% Adjustable parameters
isbackprop = 1;         % Toggle backprojection
isdiffuse  = 0;         % Toggle diffuse reflection
K          = 0;         % Downsample data to (4 ps) * 2^K = 16 ps for K = 2
isbackproj = 1;

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
% Define NLOS blur kernel 
% Compute inverse filter of NLOS blur kernel
fpsf = fftn(psf);
if (~isbackprop)
    invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./snr);
else
    invpsf = conj(fpsf);
end


% Permute data dimensions
data = permute(rect_data,[3 2 1]);

% Define volume representing voxel distance from wall
grid_z = repmat(linspace(0,1,M)',[1 N N]);

display('Inverting...');
tic;

% Step 1: Scale radiometric component
if (isdiffuse)
    data = data.*(grid_z.^4);
else
    data = data.*(grid_z.^2);
end

% Step 2: Resample time axis and pad result
tdata = zeros(2.*M,2.*N,2.*N);
tdata(1:end./2,1:end./2,1:end./2)  = reshape(mtx*data(:,:),[M N N]);

% Step 3: Convolve with inverse filter and unpad result
tvol = ifftn(fftn(tdata).*invpsf);
tvol = tvol(1:end./2,1:end./2,1:end./2);

% Step 4: Resample depth axis and clamp results
vol  = reshape(mtxi*tvol(:,:),[M N N]);
vol  = max(real(vol),0);

        if isbackproj
           % apply laplacian of gaussian filter
           vol = filterLaplacian(vol);
           
           % normalize, apply optional threshold           
           vol = max(vol./max(vol(:)),0);
        end

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




