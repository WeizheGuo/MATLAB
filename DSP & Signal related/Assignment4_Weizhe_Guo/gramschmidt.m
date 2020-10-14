function [orthset] = gramschmidt(A)
    [r,c] = size(A);
    orthset = zeros([r,c]);
    summ = zeros([r,c]);
    for i = 1:c
        for j = 1: (i-1)
            summ(:,i) = summ(:,i) + orthset(:,j)'*A(:,i)*orthset(:,j);
        end
        orthset(:,i) = A(:,i) - summ(:,i);
        orthset(:,i) = orthset(:,i)/norm(orthset(:,i));
    end
end

        




   