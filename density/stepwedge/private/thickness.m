function x_b_matrix = thickness(x_b_info,dimensions,resize_factor)

% fill a new matrix with the breast thickness at each pixel

% define the empty matrix
x_b_matrix = zeros(dimensions(1),dimensions(2));

% sort x_b_info into ascending order by xposition - this is ok, since order is no longer important
x_b_info = sortrows(x_b_info);

% get xpositions and x_b and numpoints from x_b_info
xpositions = x_b_info(:,1)*resize_factor;
x_b = x_b_info(:,2);
numpoints = size(xpositions);

if(numpoints==0) 
    x_b_matrix=input('there were no markers, please enter a thickness value in mm');
elseif(numpoints==1)
    x_b_matrix(:,:) = x_b;
else

    for i=1:dimensions(2)
        % find the pair which is to the right of this column
        if i > xpositions(numpoints(1))
            k=numpoints(1)+1;
        else
            k=1;
            while i > xpositions(k)
                k=k+1;
            end
        end
        
        % find the interpolated x_b value
        if k==1
            % extrapolate back
            x_b_interp = x_b(k) - (xpositions(k)-i)*(x_b(k+1)-x_b(k))/(xpositions(k+1)-xpositions(k));
        elseif k < (numpoints(1))+1
            % interpolate
            x_b_interp = x_b(k-1) + (i-xpositions(k-1))*(x_b(k)-x_b(k-1))/(xpositions(k)-xpositions(k-1));
        elseif k == (numpoints(1))+1
            % extrapolate forwards
            x_b_interp = x_b(k-1) + (i-xpositions(k-1))*(x_b(k-1)-x_b(k-2))/(xpositions(k-1)-xpositions(k-2));
        end
 
        % disp([i x_b_interp]);
    
        % fill this column of the matrix with this value
        this_column = ones(dimensions(1),1)*x_b_interp;
        x_b_matrix(:,i)=this_column;
    end
    
end
%     figure(100);
%     imshow(x_b_matrix);

