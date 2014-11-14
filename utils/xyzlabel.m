function xyzlabel(xlbl, ylbl, zlbl)

if exist('xlbl','var'), xlabel(xlbl); end
if exist('ylbl','var'), ylabel(ylbl); end
if exist('zlbl','var'), zlabel(zlbl); end
