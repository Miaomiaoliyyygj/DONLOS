
function [vol,tic_x,tic_y,tic_z] = cnlos_reconstruction2(rect_data,width,snr,z_trim,z_offset,psf,mtx,mtxi,bin_resolution,sampling_coeff)



% Constants
c              = 3e8;   % Speed of light (meters per second)

% Adjustable parameters
isbackprop = 0;         % Toggle backprojection
isdiffuse  = 0;         % Toggle diffuse reflection
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
 
% Define NLOS blur kernel 

% Compute inverse filter of NLOS blur kernel
fpsf = fftn(psf);
invpsf = conj(fpsf);
mtx = full(mtx);
mtxi = full(mtxi);

% Permute data dimensions
data = permute(rect_data,[3 2 1]);

% Define volume representing voxel distance from wall
grid_z = repmat(linspace(0,1,M)',[1 N N]);

display('Inverting...');
tic;

wall_size = width*2;
 % Step 0: define virtual wavelet properties
        s_lamda_limit = wall_size/(N - 1)*1.5; % sample spacing on the wall
%         sampling_coeff =1.1; % scale the size of the virtual wavelength (usually 2, optionally 3 for noisy scenes)
        virtual_wavelength = sampling_coeff * (s_lamda_limit * 2); % virtual wavelength in units of cm
        cycles = 5; % number of wave cycles in the wavelet, typically 4-6
        % Step 1: convolve measurement volume with virtual wave
        [phasor_data_cos, phasor_data_sin] = waveconv(bin_resolution, virtual_wavelength, cycles, data);
        phasor_data_cos = single(phasor_data_cos);
        phasor_data_sin = single(phasor_data_sin);
        
% Step 2: Resample time axis and pad result
  phasor_tdata_cos = single(zeros(2.*M,2.*N,2.*N));
        phasor_tdata_sin = single(zeros(2.*M,2.*N,2.*N));
        phasor_tdata_cos(1:end./2,1:end./2,1:end./2) = reshape(mtx*phasor_data_cos(:,:),[M N N]);
        phasor_tdata_sin(1:end./2,1:end./2,1:end./2) = reshape(mtx*phasor_data_sin(:,:),[M N N]);

        % Step 3: convolve with backprojection kernel
        tvol_phasorbp_sin = ifftn(fftn(phasor_tdata_sin).*invpsf);
        tvol_phasorbp_sin = tvol_phasorbp_sin(1:end./2,1:end./2,1:end./2);
        phasor_tdata_cos = ifftn(fftn(phasor_tdata_cos).*invpsf);       
        phasor_tdata_cos = phasor_tdata_cos(1:end./2,1:end./2,1:end./2);
        
        % Step 4: compute phasor field magnitude and inverse LCT
        tvol = sqrt(tvol_phasorbp_sin.^2 + phasor_tdata_cos.^2);
        vol  = reshape(mtxi*tvol(:,:),[M N N]);
vol  = max(real(vol),0);



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




