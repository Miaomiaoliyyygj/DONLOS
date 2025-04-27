function out = filterLaplacian(vol)
    % set filter parameters
    hsize = 5;
    std1 = 0.01;

    % calculate filter weights
    lim = (hsize - 1) / 2;
    std2   = std1^2;

    [x,y,z] = ndgrid(-lim:lim,-lim:lim, -lim:lim);
    w  = exp(-(x.*x + y.*y + z.*z)/(2*std2));
    w(w < eps * max(w(:))) = 0;
    sumw = sum(w(:));
    if sumw ~= 0
    w  = w/sumw;
    end;

    w1 = w.*(x.*x + y.*y + z.*z - 3*std2)/(std2^2);
    w = w1 - sum(w1(:))/hsize^3; % make the filter sum to zero

    out = convn(vol, w, 'same');
end