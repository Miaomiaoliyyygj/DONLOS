%% loading solution
close all, load(['u',num2str(loop),'.mat']), u = reshape(u,lxs,lys,lzs,3);
%% cutting boundary
% u([1:show_cut_x,end-show_cut_x+1:end],:,:,:) = [];
% u(:,[1:show_cut_y,end-show_cut_y+1:end],:,:) = [];
% u(:,:,[1:show_cut_z,end-show_cut_z+1:end],:) = [];
%% Computing albedo
albedo = sqrt(u(:,:,:,1).^2 + u(:,:,:,2).^2 + u(:,:,:,3).^2);
%% Normalization
max_albedo = max(albedo(:));
albedo = albedo / max_albedo;
u = u / max_albedo;
%% showing x y z components
y = three_components(u);
%% Three view
three(albedo, 5, loop);
saveas(gcf,'three view.png')
threshold_rate = 0;
albedo(albedo < threshold_rate) = 0;
three(albedo,6,loop);
saveas(gcf,'three view after hard-thresholding.png')