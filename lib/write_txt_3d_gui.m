% save vectors(x,y,z) in a txt file

function write_txt_3d_gui( x, y, z, name_file )
	sz_vect = length(x); %this script works with same size vectors
	fid = fopen(name_file,'wt');	
	for( i = 1 : sz_vect )
		fprintf( fid, '%f,%f,%f\n', x(i),y(i),z(i) )	;
	end 	
	fclose(fid);
end
