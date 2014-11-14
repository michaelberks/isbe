function exportfig_im(fname, im_dims, varargin)
%
% put in a filename (minus extension) and this exports
% the .fig, .eps, .pdf and .png files for use in LaTeX...
%

% export all but png by default

%Modify the figure;
set(gcf, 'windowstyle', 'normal');

fx = 1024;
fy = 1024 * im_dims(1) / im_dims(2);
set(gcf, 'units', 'pixels', 'position', [10 10 fx fy]); 
set(gca, 'units', 'normalized', 'position', [0 0 1 1]);

opts = {'fig','eps','pdf','png'};
if length(varargin)>0, opts = varargin{1}; end

[p,n,e] = fileparts(fname);
if isempty(p)
	p = './';
end

if ~exist(p,'dir')
	mkdir(p);
end

if (strcmp(e,'.fig') | strcmp(e,'.eps') | ...
		strcmp(e,'.pdf') | strcmp(e,'.png'))
	fname = fullfile(p,n);
end

fprintf('Exporting...');

errmsg = '';

% save figure
if any(strcmp(opts(:),'fig'))
	fprintf('fig...');
	saveas(gcf,[fname '.fig']);
end

% export eps
if any(strcmp(opts(:),'eps'))
	fprintf('eps...');
% 	print('-depsc2','-loose',[fname '.eps']);
	print('-depsc2', [fname '.eps']);
end

% export png
if any(strcmp(opts(:),'png'))
	fprintf('png...');
	try
		eval(['!convert -density 150 "',fname,'.eps" "',fname,'.png"']);
	catch
		print('-dpng', [fname '.png']); % ugly output - ImageMagick is better
	end
end

% convert to pdf using epstopdf
if any(strcmp(opts(:),'pdf'))
	fprintf('pdf...');
	try
		if (isunix)
			eval(['!epstopdf ' fname '.eps']);
		else
            fontpath = 'C:\Windows\Fonts\';
%             epstopdf_cmd = ['epstopdf -sFONTPATH="', fontpath, '" '];
            epstopdf_cmd = ['epstopdf '];
			if ~isempty(p)
				olddir	= cd;
				cd(p);
					eval(['!', epstopdf_cmd, ' ', n, '.eps']);
				cd(olddir);	
            else
    			eval(['!', epstopdf_cmd, ' ', n, '.eps']);
			end
		end
	catch
		errmsg = [errmsg,'Cannot export PDF - install epstopdf\n'];
	end
end
set(gcf, 'windowstyle', 'docked');
fprintf('done\n');
fprintf('%s',errmsg);
