function box = drawbox(boxsize,dim)
    switch dim
        case 3
            [xx,yy,zz] = meshgrid(-5:5,-5:5,-5:5);
            box = exp(-(xx.^2+yy.^2+zz.^2)/(boxsize^2));
            box = box/sum(box(:));
        case 2
            [xx,yy] = meshgrid(-5:5,-5:5);
            box = exp(-(xx.^2+yy.^2)/(boxsize^2));
            box = box/sum(box(:));
            box = reshape(box,[1,11,11]);
            
    end
end