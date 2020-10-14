function matrix = InsertZeros(x)
    matrix = zeros(size(x,1)*2, size(x,2)*2);
    for row = 1 : size(x,1)
        for col = 1 : size(x,2)
            matrix(row*2, col*2) = x(row, col);
        end
    end
    matrix = matrix(2:end,2:end);
end
