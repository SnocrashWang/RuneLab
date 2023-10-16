function [y] = meanFilter(x, size)
    b = ones(1,size)./size;
    y = conv(x,b);
    y = y(size : length(y)-size+1);
    y = [x(1:fix(size/2)), y, x(length(x)-fix(size/2)+1:length(x))];
end
