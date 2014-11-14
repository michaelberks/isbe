function window_motion_test
%
figure('WindowButtonDownFcn',@wbdcb)
ah = axes('DrawMode','fast');
axis([1 10 1 10]); axis equal; hold all;
title('Click and drag')
xinit = [];
yinit = [];
    function wbdcb(src,evnt)
         if strcmp(get(src,'SelectionType'),'normal')
            set(src,'pointer','circle')
            cp = get(ah,'CurrentPoint');
            xinit = cp(1,1);yinit = cp(1,2);
    %         hl = line('XData',xinit,'YData',yinit,...
    %         'Marker','p','color','b');
            hl = plot(xinit, yinit, 'b-x', 'XDataSource','xinit','YDataSource','yinit');
            set(src,'WindowButtonMotionFcn',@wbmcb)
            set(src,'WindowButtonUpFcn',@wbucb)
         end
    end
    function wbmcb(src,evnt)
       cp = get(ah,'CurrentPoint');
       xinit = [xinit,cp(1,1)];
       yinit = [yinit,cp(1,2)];
       %set(hl,'XData',xdat,'YData',ydat);
       refreshdata(ah, 'caller');
       %drawnow
    end
   
    function wbucb(src,evnt)
       if strcmp(get(src,'SelectionType'),'alt')
          set(src,'Pointer','arrow')
          set(src,'WindowButtonMotionFcn','')
          set(src,'WindowButtonUpFcn','')
       else
           plot(xinit, yinit, 'r');
          return
       end
    end
end