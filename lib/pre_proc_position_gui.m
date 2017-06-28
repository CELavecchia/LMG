%clc; clear all; close all;
function [Lmes,VTstr] = pre_proc_position_bdynFU( Lmes, VTstr, nbodies)
%function [L2,L3,I1,I2] = pre_proc_position_bdynFU(L1,L2,L3,I1,I2)


%function [L2,L3,I1] = pre_proc_position(L1,L2,L3,I1)
if(nbodies == 2)
    
    %% distance between L1 and IVD1
    L1dist = Lmes(1).VT(Lmes(1).VT(:,2)>-20,:);

    d1 = min(L1dist(:,3))-max(VTstr.VT(:,3))+3;
    t1 = [ 0 0 d1];
    VTstr.VT(:,3) = VTstr.VT(:,3)+d1;
    %plot_matrix(I1,'.g')

    %% distance between IVD1 and L2
    L2dist = Lmes(2).VT(Lmes(2).VT(:,2)>-20,:);
    %plot_matrix(L2dist,'.y');

    d2 = min(VTstr.VT(:,3)) - max(L2dist(:,3))+3;

    t2 = [ 0 0 d2];%+7];
    Lmes(2).VT(:,3) = Lmes(2).VT(:,3)+t2(:,3);
%plot_matrix(L2,'.r');

    %% distance between IVD1 and L2
    L2dist = Lmes(2).VT(Lmes(2).VT(:,2)>-20,:);
    %plot_matrix(L2dist,'.y');

    d2 = min(VTstr.VT(:,3)) - max(L2dist(:,3))+3;

    t2 = [ 0 0 d2];%+7];
    Lmes(2).VT(:,3) = Lmes(2).VT(:,3)+t2(:,3);
    %plot_matrix(L2,'.r');
elseif nbodies ==3

    %% distance between L2 and IVD2
    d3 = min(L2dist(:,3)) - max(VTstr.VT2(:,3))+3;

    t3 = [ 0 0 d3];%+1.5];
    VTstr.VT2(:,3) = VTstr.VT2(:,3)+t3(3);
    
    
    %% distance between IVD2 and L3
    L3dist = Lmes(3).VT(Lmes(3).VT(:,2)>-20,:);
    %plot_matrix(L2dist,'.y');
    d4 = min(VTstr.VT2(:,3)) - max(L3dist(:,3))+3;

    t4 = [ 0 1 d4];%+6];%1.5];
    Lmes(3).VT = Lmes(3).VT+t4;
    %plot_matrix(L3,'.r');
    %distance between L2 and IVD2
end
end




