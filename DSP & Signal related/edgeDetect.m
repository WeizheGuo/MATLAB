function bool = edgeDetect(g,T, order)
    hy = 1/8*[-1,-2,-1;0,0,0;1,2,1];
    hx = hy';
    xconv = conv2(g, hx,'same');
    yconv = conv2(g, hy,'same');
    xconv = xconv(1:size(g,1),1:size(g,2));
    yconv = yconv(1:size(g,1),1:size(g,2));
    if order ==2
        result = sqrt(xconv.^2 + yconv.^2);
    elseif order==1
        result = abs(xconv) + abs(yconv);
    end
    T_m = T*ones(size(g,1), size(g,2));
    bool = (result>T_m);
end