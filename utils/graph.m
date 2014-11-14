function H = graph(varargin)

H = gcf;
if length(varargin)>0, H = figure(varargin{1}); end

% if (isunix),	figpos = [960 720 640 480];
% else,					figpos = [920 600 480 360];
% end

% set box to on for all axes
h_child = get(H,'Children');
for i = 1:length(h_child)
	if strcmp(get(h_child(i),'type'),'axes')
		set(h_child(i),'Box','on'); 
	end
end

% set(H,'PaperPosition',[0 0 12 9],...
% 			'Renderer','Painters',...
% 			'Position',figpos,...
% 			'DoubleBuffer','on');
		
set(H,'PaperPosition',[0 0 9 9],...
			'Renderer','Painters',...
			'DoubleBuffer','on');

