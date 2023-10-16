function [rmse] = RMSE(x)
    rmse = sqrt(sum((x).^2) ./ length(x));
end
