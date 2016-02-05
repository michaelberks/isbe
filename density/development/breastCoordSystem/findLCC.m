function LCC = findLCC(binaryIm, connectivity)
  % Find the largest connected component from a binary classification

  [ACC, nbr] = bwlabeln(binaryIm, connectivity); % All Connected Components
  if ~nbr
     disp('Found no components')
     LCC = 0*ACC;
     return
  end

  % Find the largest connected component
  compCount = zeros(1,nbr);
  for c = 1:nbr
     compCount(c) = sum(ACC(:)==c);
  end

  [count,comp] = max(compCount);
  LCC = ACC==comp; % Binary for largest connected component

%   if ~isempty(globalOptions)
%      ratio = Get(globalOptions,'LCCratio');
%      [counts, comps] = sort(compCount);
%      for i = (nbr-1):-1:1 % last is largest
%         c = counts(i);
%         if c > count*ratio
%            fprintf('Adding extra component, voxel count %d -> %d \n',count,count+c);
%            comp = comps(i);
%            count = count + c;
%            LCC = LCC | ACC==comp;
%         end
%      end
%   end
  
