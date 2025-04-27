function y = three( albedo, fig_num, j)

min_albedo = 0; 
max_albedo = max(albedo(:));

figure(fig_num)
% front view
    subplot(1,3,1)
    temp_img = rot90(squeeze(max(albedo,[],1)));
    imagesc(temp_img(:,end:-1:1));
    xlabel('y');
    ylabel('z');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title(['Front view: iter ',num2str(j)])
    axis off

% top view  
    subplot(1,3,2)
    temp_img = squeeze(max(albedo,[],3));
    imagesc(temp_img(end:-1:1,end:-1:1));
    xlabel('y');
    ylabel('x');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title(['Top view: iter ',num2str(j)])
    axis off
    
% side view
    subplot(1,3,3)
    temp_img = rot90(squeeze(max(albedo,[],2)));
    imagesc(temp_img);
    xlabel('x');
    ylabel('z');
    colormap('gray');
    caxis([min_albedo, max_albedo])
    axis equal;
    title(['Side view: iter ',num2str(j)])
    axis off

   
    y = [];
end