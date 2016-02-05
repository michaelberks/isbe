function showPoints(objCS)

    hold on;

    plot(objCS.xA,objCS.yA,'ro')
    hold on;
    plot(objCS.xA,objCS.yA,'kx')
    plot(objCS.xB,objCS.yB,'ro')

    plot(objCS.xB,objCS.yB,'rx')
    plot(objCS.xC,objCS.yC,'ko')
    plot(objCS.xC,objCS.yC,'rx')

    if ~isempty(objCS.yM)
        plot(objCS.xM,objCS.yM,'gx')
        plot(objCS.xM,objCS.yM,'go')    
    end
    
    h = text(objCS.xB,objCS.yB,'B','fontsize',10,'color','green');
    text(objCS.xA,objCS.yA,'A','fontsize',10,'color','green');
    text(objCS.xC,objCS.yC,'C','fontsize',10,'color','green');

    xs = [objCS.xA objCS.xB objCS.xC];
    ys = [objCS.yA objCS.yB objCS.yC];

%     xlim([max(-160,min(xs)-160) max(xs)+100])
%     ylim([max(-160,min(ys)-160) max(ys)+100])
    
    %title('Breast Coord System')
    hold off;
end

