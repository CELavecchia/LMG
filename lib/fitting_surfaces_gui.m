% function to fit the upper and inferior surfaces of the vertebrae to
% replicate the same geometry for the disc
function interpFunc = fitting_surfaces_gui(EP);

    for(j=1:5)
       interpFunc(j).sup =  scatteredInterpolant(EP(j).sup(:,[1 2]),EP(j).sup(:,3),'natural');
       interpFunc(j).inf =  scatteredInterpolant(EP(j).inf(:,[1 2]),EP(j).inf(:,3),'natural');
    end

end