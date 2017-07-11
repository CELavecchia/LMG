%{
 function to save points on the superior and inferior surfaces. 
%}
function save_coord_gui(L_no,cloudpathName);

ind = 1 ;
    for(j = 1:5)
        fprintf('saving the coordinates for the vertebrae L%d \n',j)
            
        outputName=fullfile(cloudpathName,sprintf('L_body_%d.txt',j));
        write_txt_3d_gui( L_no(j).body(:,1), L_no(j).body(:,2), L_no(j).body(:,3), outputName );
        
        outputName=fullfile(cloudpathName,sprintf('L_lam_%d.txt',j));
        write_txt_3d_gui( L_no(j).lam(:,1), L_no(j).lam(:,2), L_no(j).lam(:,3), outputName );
        
        
        outputName=fullfile(cloudpathName,sprintf('L_ped_%d.txt',j));
        write_txt_3d_gui( L_no(j).ped(:,1), L_no(j).ped(:,2), L_no(j).ped(:,3), outputName );

        outputName=fullfile(cloudpathName,sprintf('L_proc_grid_%d.txt',j));
        write_txt_3d_gui( L_no(j).proc(:,1), L_no(j).proc(:,2), L_no(j).proc(:,3), outputName );
        
        ind = ind+3;
    end


end

