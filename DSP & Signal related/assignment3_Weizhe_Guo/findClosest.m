function[val, ind] = findClosest(x,desiredValue)
    y = x(:);
    y1 = abs(y-desiredValue);
    val = min(y1)+desiredValue;
    ind = find(min(y1)==y1);
end

