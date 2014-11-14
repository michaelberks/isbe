function n = make_neighbourhood2(scale)
% MAKE_NEIGHBOURHOOD - 
%   
  
  scale = scale + 2;
  n = zeros(scale);
  
  mid = ceil(scale / 2);
  mid = [mid, mid];
  for i=1:scale,
    for j=1:scale,
      d = [i, j] - mid;
      d = sqrt(sum(d.^2));
      if(d <= (scale / 2))
	n(i, j) = 1;
      end
    end
  end
  