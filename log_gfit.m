
%cd('/home/sgrlab/Documents/CJK/PEGDA_HGS/')
cd('/Users/cjk/Box/CJ_Kelsheimer/PEGDA_HGS')

fname = 'pegda_globalfit_tracker.txt';
fid = fopen(fname,'a');

fprintf(fid,'\n');
fprintf(fid,'----\n');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'\n');

fprintf(fid,'%s %s %s %s %s %s %s %s %s\n',options.pnames{:});
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'Original parameters\n');
fprintf(fid,'%f ',p0);
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'Lower bounds\n');
fprintf(fid,'%f ',lb);
fprintf(fid,'\n');

fprintf(fid,'Upper bounds\n');
fprintf(fid,'%f ',ub);
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'Output (fit) parameters\n');
fprintf(fid,'%f ',pfit);
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'Minimized Square Error = %f',err);
fprintf(fid,'\n');


fclose(fid);

