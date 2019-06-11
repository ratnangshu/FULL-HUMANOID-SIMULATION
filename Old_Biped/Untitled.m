Rh=zeros(4 , 6);
lh=zeros(4 , 6);
rha =zeros(4,1);
lha =zeros(4,1);
waist=zeros(4,4);
rleg=zeros(4,8);
lleg=zeros(4,8);
rlegz=zeros(200,8);
llegz=zeros(200,8);
l0=zeros(7,1);
r0=zeros(7,1);
rfoot=zeros(4,8);
lfoot=zeros(4,8);
xcom=zeros(4,1);
z=30*cos(pi*18/180);
h=z;
ll = 15;
y=0;
t=0;
t0 = [1 0 0 0;0 1 0 0; 0 0 1 0;0 0 0 1];
rf = [1 0 0 0;0 0 1 0; 0 -1 0 0;0 0 0 1];
wa1=0.1;
wa2=0.01;
wa = 0.05;
sl=16;
z__1=zeros(4,1);
z_1=zeros(4,1);
la=zeros(7,200);
ra=zeros(7,200);
waist_1=zeros(80,1);
w=0;



    %###############################################################
    %######    right leg as swing leg                          #####
    %###############################################################
    %{
    for t=1:40
        if (t>-1) && (t<21)
            x=6;
            xa=-5;
            y= (3*t)/20;
            ha=0;
        else
            x=6;
            xa=-5;
            y= 3;
            ha=(2.25/20)*(t-20);
        end
  
        l0(6)=atan((-y)/h);
        l0(2)=-l0(6);
        [l0(3), l0(4), l0(5)] = calc_angles(-5,0,x,h,ll);
        l0(5)=-l0(5);
        l0(4)=-l0(4);
        
        
        r0(6) = l0(6);
        r0(2) = l0(2);
        [r0(3), r0(4), r0(5)] = calc_angles(xa,ha,x,h,ll);
        r0(4)=-r0(4);
        r0(5)=-r0(5);
        
        
        rf(:,4)=lleg(:,7);
        [lleg, lfoot, xcom] = calc_pos(l0,rf,1);
        t0(:,4)=xcom;
        [rleg, rfoot] = calc_rpos(r0,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;
        plot_final(rleg,lleg,waist,zmp);
        la(:,t)=l0;
        ra(:,t)=r0;
        pause(.01);
    end
  
    for t=60:100
        x=t/10;
        if (t>59) && (t<81)
            y = 3;
            xa = (20/4)*(x-4) - 5;
            ha = 0.00237671*(xa)^3 -0.0582294*(xa)^2 + 0.16637*(xa) + 2.58467;
        else
            y = -3*(t-80)/20 + 3;
            xa = 15;
            ha = 0;
        end


        l0(6)=atan((-y)/h);
        l0(2)=-l0(6);
        [l0(3), l0(4), l0(5)] = calc_angles(-xa,ha,x,h,ll);
        l0(5)=-l0(5);
        l0(4)=-l0(4);

        r0(6) = l0(6);
        r0(2) = l0(2);
        [r0(3), r0(4), r0(5)] = calc_angles(-5,0,x,h,ll);
        r0(4)=-r0(4);
        r0(5)=-r0(5);

        swap = r0;
        r0 = l0;
        l0 = swap;

        la(:,t-20)=l0;
        ra(:,t-20)=r0;

        %position calculation ---------------------------------------------
        rf(:,4)=lleg(:,7);
        [lleg, lfoot, xcom] = calc_pos(l0,rf,1);
        t0(:,4)=xcom;
        [rleg, rfoot] = calc_rpos(r0,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;
        plot_final(rleg,lleg,waist,zmp);
        pause(0.01);
        
    end
    %}
    for t=1:40
        %t=0;

        if (t>-1) && (t<11)
            x = t/5;
            y =-(t)/2;
            xa = -3;
            ha = 0;
            hb =h;
        elseif (t>10) && (t<31)
            x = 1 + t/10;
            y = -5;
            xa = (3/5)*t - 9;
            ha = -0.005*(xa+20)*(xa+3)*(xa-9);
            hb= -0.04*(x-4)*(x-2)*(x-1.2) + h;
        else
            x = t/5 -2;
            y = (t-30)/2-5;
            xa = 9;
            ha = 0;
            hb=h;
        end
        
        l0(6)=atan((-y)/h);
        l0(2)=-l0(6);
        [l0(3), l0(4), l0(5)] = calc_angles(-xa,ha,x,hb,ll);
        l0(5)=-l0(5);
        l0(4)=-l0(4);



        r0(6) = l0(6);
        r0(2) = l0(2);
        [r0(3), r0(4), r0(5)] = calc_angles(-3,0,x,hb,ll);
        r0(4)=-r0(4);
        r0(5)=-r0(5);

        la(:,t)=l0;
        ra(:,t)=r0;

        %position calculation ---------------------------------------------
        rf(:,4)=rleg(:,7);
        [rleg, rfoot, xcom] = calc_pos(r0,rf,-1);
        t0(:,4)=xcom;
        [lleg, lfoot] = calc_lpos(l0,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;
        zmp(1)=zmp(1)-2;
        waist_1(w+1)=waist(3,2);
        w=w+1;
        [z__1, z_1] = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1);
        pause(0.1); %for this just add pause(wa) and then comment it.
    end

    %###############################################################
    %######    left leg as swing leg                          #####
    %###############################################################

    for t=1:40
        %t=0;
       
        if (t>-1) && (t<11)
            x = t/5;
            y =(t)/2;
            xa = -3;
            ha = 0;
            hb =h;
        elseif (t>10) && (t<31)
            x = 1 + t/10;
            y = 5;
            xa = (3/5)*t - 9;
            ha = -0.005*(xa+20)*(xa+3)*(xa-9);
            hb= -0.04*(x-4)*(x-2)*(x-1.2) + h;
        else
            x = t/5 -2;
            y = -(t-30)/2+5;
            xa = 9;
            ha = 0;
            hb=h;
        end


        l0(6)=atan((-y)/h);
        l0(2)=-l0(6);
        [l0(3), l0(4), l0(5)] = calc_angles(-xa,ha,x,hb,ll);
        l0(5)=-l0(5);
        l0(4)=-l0(4);

        r0(6) = l0(6);
        r0(2) = l0(2);
        [r0(3), r0(4), r0(5)] = calc_angles(-3,0,x,hb,ll);
        r0(4)=-r0(4);
        r0(5)=-r0(5);

        swap = r0;
        r0 = l0;
        l0 = swap;

        la(:,t+40)=l0;
        ra(:,t+40)=r0;

        %position calculation ---------------------------------------------
        rf(:,4)=lleg(:,7);
        [lleg, lfoot, xcom] = calc_pos(l0,rf,1);
        t0(:,4)=xcom;
        [rleg, rfoot] = calc_rpos(r0,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;
        zmp(1)=zmp(1)-2;
        waist_1(w+1)=waist(3,2);
        w=w+1;
        [z__1, z_1] = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1);
       
        pause(0.1);
    end
for i=1:80
    %la(5,i)=5.9810-la(5,i);
    %ra(2,i)=5.1695-ra(2,i);
    %la(3,i)=3.2290+la(3,i);
    %ra(5,i)= -1.2210+ra(5,i);
    %ra(4,i)=1.3805-ra(4,i);
    %la(4,i)=2.9651+la(4,i);
    %la(2,i)=4.8059-la(2,i);
    %la(5,i)=4.1110+la(5,i);
    fprintf('[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f],',-la(3,i), -ra(2,i), -la(5,i), -la(2,i), ra(5,i), -ra(4,i), -la(4,i), -ra(2,i), -la(2,i), ra(3,i));
end