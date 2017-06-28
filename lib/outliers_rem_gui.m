% removing outliers coming from the fitting 
function L_body_c = outliers_rem(L_body)
%remove outliers from vert bodies
%count = 1;
count2 = 1;
for( j = 1:5 )
    count = 1;
    for( i = 1:length(L_body) )
         if( L_body(i,count2+1)<= -20.00 && L_body(i,count2+2)<= -15.0)  %modified on the 12/01/17
            exclude = L_body(i,:);
         else
            L_body_c(count,count2:count2+2) = L_body(i,count2:count2+2);
            count = count+1;
         end
    end
    count2 = j*3+1;
end

end