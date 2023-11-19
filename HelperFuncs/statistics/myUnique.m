function y = myUnique(x)
  y = unique(x);
%   y(isnan(y(1:end-1))) = [];
  y(isnan(y)) = [];
end