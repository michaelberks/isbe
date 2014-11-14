% fix random seed
rand('twister',314159);

n_by_3 = 1600;
n = n_by_3 * 3;
s = 10^ceil(log10(n));
inds = zeros(n,3);

outpath = 'S:\projects\nailfold\image_lists\randomized';

for suffix = {'c','d'}
    row1 = 1;
    for i = 1:3
        c1 = randperm(s)-1;
        c1 = i*s + c1(1:n_by_3)';

        c2 = randperm(s)-1;
        c2 = (i+3)*s + c2(1:n_by_3)';

        c3 = randperm(s)-1;
        c3 = (i+6)*s + c3(1:n_by_3)';

        inds(row1:row1+n_by_3-1,:) = [c1,c2,c3];
        row1 = row1 + n_by_3;
    end

    % shuffle columns
    for i = 1:n
        switch ceil(rand*3)
            case 1, inds(i,[1,2,3]) = inds(i,[1,2,3]);
            case 2, inds(i,[1,2,3]) = inds(i,[2,3,1]);
            case 3, inds(i,[1,2,3]) = inds(i,[3,1,2]);
        end
    end

    % shuffle rows
    inds = inds(randperm(n),:);
    
    u = unique(inds);
    [numel(u) numel(inds)]
    
    for j = 1:3
        filename = sprintf('list%i_%c.txt', j, suffix{1});
        fid = fopen([outpath,'/',filename], 'w');
        for i = 1:n
            fprintf(fid, '%d%c\n', inds(i,j), suffix{1});
        end
        fclose(fid);
    end
end
