function y = three_components(u)
[m,n,o,~] = size(u);
pos = NaN * zeros(n,o);
com_z = zeros(n,o);
com_y = zeros(n,o);
com_x = zeros(n,o);
for j = 1:n
    for k = 1:o
        temp = u(:,j,k,1); temp = temp(:); 
        if sum(abs(temp)) == 0, continue, end
        f = find(temp == min(temp)); 
        pos(j,k) = f(1);
        com_z(j,k) = -u(pos(j,k),j,k,1);
        com_y(j,k) = -u(pos(j,k),j,k,3);
        com_x(j,k) = -u(pos(j,k),j,k,2);
    end
end
com_x = rot90(com_x); com_x = com_x(:,end:-1:1);
com_y = rot90(com_y); com_y = com_y(:,end:-1:1);
com_z = rot90(com_z); com_z = com_z(:,end:-1:1);

% x component
temp_img = com_x; min_temp_img = min(temp_img(:)); max_temp_img = max(temp_img(:)); 
figure('Position',[50,50,320,400]);
axes('Position',[0.02 0.212 0.96 0.768]);
imagesc(temp_img), colormap('hot'), axis equal, axis off

c = colorbar;
c.Location = 'south';
c.Position = [0.02 0.11 0.96 0.08];
c.Limits = [min_temp_img,max_temp_img];
c.TickLength = 0;
c.Ticks = [min_temp_img + 0.09 * (max_temp_img - min_temp_img ) ,max_temp_img -  0.09 * (max_temp_img - min_temp_img )];
c.TickLabels = {num2str(min_temp_img,'%.2f'),num2str(max_temp_img,'%.2f')};
c.FontName = 'Times New Roman';
c.FontSize = 23;
c.AxisLocation = 'out';
c.Box = 'on';

saveas(gcf,'x-component.png')

% y component
temp_img = com_y; min_temp_img = min(temp_img(:)); max_temp_img = max(temp_img(:)); 
figure('Position',[50,50,320,400]);
axes('Position',[0.02 0.212 0.96 0.768]);
imagesc(temp_img), colormap('hot'), axis equal, axis off

c = colorbar;
c.Location = 'south';
c.Position = [0.02 0.11 0.96 0.08];
c.Limits = [min_temp_img,max_temp_img];
c.TickLength = 0;
c.Ticks = [min_temp_img + 0.09 * (max_temp_img - min_temp_img ) ,max_temp_img -  0.09 * (max_temp_img - min_temp_img )];
c.TickLabels = {num2str(min_temp_img,'%.2f'),num2str(max_temp_img,'%.2f')};
c.FontName = 'Times New Roman';
c.FontSize = 23;
c.AxisLocation = 'out';
c.Box = 'on';

saveas(gcf,'y-component.png')

% z component
temp_img = com_z; min_temp_img = min(temp_img(:)); max_temp_img = max(temp_img(:)); 
figure('Position',[50,50,320,400]);
axes('Position',[0.02 0.212 0.96 0.768]);
imagesc(temp_img), colormap('hot'), axis equal, axis off

c = colorbar;
c.Location = 'south';
c.Position = [0.02 0.11 0.96 0.08];
c.Limits = [min_temp_img,max_temp_img];
c.TickLength = 0;
c.Ticks = [min_temp_img + 0.09 * (max_temp_img - min_temp_img ) ,max_temp_img -  0.09 * (max_temp_img - min_temp_img )];
c.TickLabels = {num2str(min_temp_img,'%.2f'),num2str(max_temp_img,'%.2f')};
c.FontName = 'Times New Roman';
c.FontSize = 23;
c.AxisLocation = 'out';
c.Box = 'on';

saveas(gcf,'z-component.png')

%% albedo
temp_albedo = sqrt(com_x .^2 + com_y.^2 + com_z.^2);

temp_img = temp_albedo; min_temp_img = min(temp_img(:)); max_temp_img = max(temp_img(:)); 
figure('Position',[50,50,320,400]);
axes('Position',[0.02 0.212 0.96 0.768]);
imagesc(temp_img), colormap('hot'), axis equal, axis off

c = colorbar;
c.Location = 'south';
c.Position = [0.02 0.11 0.96 0.08];
c.Limits = [min_temp_img,max_temp_img];
c.TickLength = 0;
c.Ticks = [min_temp_img + 0.09 * (max_temp_img - min_temp_img ) ,max_temp_img -  0.09 * (max_temp_img - min_temp_img )];
c.TickLabels = {num2str(min_temp_img,'%.2f'),num2str(max_temp_img,'%.2f')};
c.FontName = 'Times New Roman';
c.FontSize = 23;
c.AxisLocation = 'out';
c.Box = 'on';

saveas(gcf,'albedo.png')

y =[];
end