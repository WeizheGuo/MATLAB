function v = orthoProj(w,ONs)
    ONS = gramschmidt(ONs);
    [~,c] = size(ONS);
    [r,c2] = size(w);
    v = zeros([r,c2]);
    for i = 1:c2
        for j = 1:c
        v(:,i) =  v(:,i) + (ONS(:,j)'*w(:,i))*ONS(:,j);
        end
    end
end
