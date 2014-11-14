% script to convert first attempt at a markup to newer format

datadir = 'U:\projects\nailfold\images\marked\';
d = dir([datadir,'*_markup.txt']);

for i = 1%:d
	filename_in = [datadir,d(i).name];
	vessels = read_vessels_from(filename_in);
	if isempty(vessels)
		
	
	filename_out = strrep(filename_in,'markup.txt','markup.000.txt');
	movefile(filename_in, filename_out);
	
	fid = fopen(filename_in,'w');
		fprintf(fid,'ncm_annotation: {\n');
		fprintf(fid,'  version: 1\n');
		fprintf(fid,'  observer: -\n');
		fprintf(fid,'  timestamp: 0\n');
		fprintf(fid,'  time_taken: -1\n');
		fprintf(fid,'  vessels: {\n');
		
		for v = 1:length(vessels)
			fprintf(fid,'    ncm_vessel: {\n');
			fprintf(fid,'      version: 1\n');
			fprintf(fid,'      points: {\n');
			fprintf(fid,'        // (x, y, width)\n');
			for p = 1:size(vessels{v},1)
				fprintf(fid,'        %d %d 0\n',vessels{v}(p,1),vessels{v}(p,2));
			end
			fprintf(fid,'      } // points\n');
			fprintf(fid,'      apex_indices: {\n');
			fprintf(fid,'      } // apex_indices\n');
			fprintf(fid,'    } // ncm_vessel\n');
		end
		
		fprintf(fid,'  }\n');
		fprintf(fid,'}\n');
	fclose(fid);	
end