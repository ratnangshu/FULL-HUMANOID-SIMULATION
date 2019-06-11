 function [lleg, lfoot] = calc_lpos(l0,t0)
%left leg ----------------------------------------------
        lt11 = ([0 -1 0 0;-1 0 0 5;0 0 -1 0;0 0 0 1])*Rz(l0(1));
        lt1 = ([0 -1 0 0;0 0 -1 0;1 0 0 0;0 0 0 1])*Rz(l0(2));
        lt2 = ([1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1])*Rz(l0(3));
        lt3 = ([1 0 0 15;0 1 0 0;0 0 1 0;0 0 0 1])*Rz(l0(4));
        lt4 = ([1 0 0 15;0 1 0 0;0 0 1 0;0 0 0 1])*Rz(l0(5));
        lt5 = ([1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 1])*Rz(l0(6));
        lt6 = ([0 1 0 0;0 0 1 0;1 0 0 5;0 0 0 1])*Rz(l0(7));
        %--------------------------------- %left leg ---------------------------------------------
        lleg(:,1) = t0*lt11*([0;0;0;1]);
        lleg(:,2) = t0*lt11*lt1*([0;0;0;1]);
        lleg(:,3) = t0*lt11*lt1*lt2*([15;0;0;1]);
        lleg(:,4) = t0*lt11*lt1*lt2*lt3*([15;0;0;1]);
        lleg(:,5) = t0*lt11*lt1*lt2*lt3*lt4*([0;0;0;1]);
        lleg(:,7) = t0*lt11*lt1*lt2*lt3*lt4*lt5*([0;0;5;1]);
        lleg(:,8) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([5;0;0;1]);
        lleg(:,6) = t0*lt11*lt1*lt2*lt3*lt4*lt5*([0;0;-5;1]);
        %--------------------------------------------------------------------------------------
        %left foot-------------------------------------------------
        lfoot(:,1) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([0;0;-4.5;1]);
        lfoot(:,2) = t0*lt11*lt1*lt2*lt3*lt4*lt5*([0;-4.5;-5;1]);
        lfoot(:,3) = t0*lt11*lt1*lt2*lt3*lt4*lt5*([0;4.5;-5;1]);
        lfoot(:,4) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([0;0;4.5;1]);
        lfoot(:,5) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([5;0;4.5;1]);
        lfoot(:,6) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([5;0;-4.5;1]);
        lfoot(:,7) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([0;0;-4.5;1]);
        lfoot(:,8) = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*([0;0;4.5;1]);
 end
