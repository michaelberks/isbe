function plot_green(CPS)

    [xx, yy] = meshgrid(-2:.05:2, -2:.05:2);
    resq = xx.^2 + yy.^2;
    resq(resq == 0) = 1;
    if CPS
        zz = -resq.*log(resq);
    else
        A = 1 ./ resq;
        zz = -resq.*(A - 1 - log(A));
        greens(1, 3)
    end

    figure, surf(xx, yy, zz);
    hold on;
    
    %%%%%%%%%%%%%%%%
    function z = greens(x, y)
        z = x + y;
    end
end