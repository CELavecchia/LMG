function plot_matrix(matrix,c); 
    
   plot3(matrix(:,1),matrix(:,2),matrix(:,3),c),hold on,grid on,xlabel('x'),ylabel('y'),zlabel('z'); 

end