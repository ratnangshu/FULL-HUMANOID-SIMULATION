% workout on yourself for in and outputs
% This helps in plotting of the data points and returns some past ZMP position for calc of acceleration and velocity.

function [z2, z1] = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1)
        plot3(rleg(1,:),rleg(2,:),rleg(3,:),'k -o',lleg(1,:),lleg(2,:),lleg(3,:),'k -o',waist(1,:),waist(2,:),waist(3,:),'k -o',zmp(1),zmp(2),zmp(3),'b -o',rfoot(1,:),rfoot(2,:),rfoot(3,:),'k -o',lfoot(1,:),lfoot(2,:),lfoot(3,:),'k -o','LineWidth',1.5);
        
               
        %az = azimuthal angle
        %el = elevation
        
        %az = 0;
        %el = 0;
        %view(az, el);
        
        %xlim([-30 30])
        %ylim([-15 15])
        %zlim([0 60])
                
        %view([0 0]);
        z=30*cos(pi*15/180);
        d=z/(980*0.000001);
        a=-d;
        b=2*d+1;
        c=-d;

        p = a.*z__1 + b.*z_1 + c.*zmp;
        %plot(z_1(1),z_1(2),'b--o');
        %plot(zmp(1),zmp(2),'b--o');
        %if (p(2)<10) && (p(2)>-10)
         %   plot(p(1),p(2),'r--o');
        %end
        grid on;
        box on;
        axis equal;
        %axis square;
        %axis([-15 150 -15 15 -15 50]);
        hold off;
        z2 = z_1;
        z1 = zmp;
end
