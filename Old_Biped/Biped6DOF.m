%##### Variables Declarartion of appropriate size ##############

Rh=zeros(4,6);
lh=zeros(4,6);
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
wa1=0.001;
wa2=0.001;
wa = 0.01;
sl=16;
z__1=zeros(4,1);
z_1=zeros(4,1);
la=zeros(7,200);
ra=zeros(7,200);
count = 1;
waist_val = zeros(200,1);



    %###############################################################
    %######    right leg as swing leg                          #####
    %###############################################################

    for t=1:100
        %t=0;

        x=(t)/10;

        if (t>-1) && (t<41) % SSP1
            y=-(3*t)/40;
            xa = -5;
            ha = 0;
        elseif (t>40) && (t<81) % DSP1
            y = -3;
            xa = (20/4)*(x-4) - 5;
            ha = 0.00237671*(xa)^3 -0.0582294*(xa)^2 + 0.16637*(xa) + 2.58467;
        else % SSP2
            y = 3*(t-80)/20 - 3;
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

        la(:,t)=180.*l0/pi;
        ra(:,t)=180.*r0/pi;

        %position calculation ---------------------------------------------
        rf(:,4)=rleg(:,7);
        [rleg, rfoot, xcom] = calc_pos(r0,rf,-1);
        t0(:,4)=xcom;
        [lleg, lfoot] = calc_lpos(l0,t0);
        lh=cal_lh(lha,t0);
        Rh=cal_rh(rha,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;

        waist_val(count,1)=waist(3,1);
        count = count + 1;
        
        [z__1, z_1] = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1);
        rlegz(t,:) = rleg(3,:);
        llegz(t,:) = lleg(3,:);
        pause(wa); %for this just add pause(wa) and then comment it.
    end

    %###############################################################
    %######    left leg as swing leg                          #####
    %###############################################################

    for t=1:100
        %t=0;
        x=(t)/10;

        if (t>-1) && (t<41)
            y=(3*t)/40;
            xa = -5;
            ha = 0;
        elseif (t>40) && (t<81)
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

        la(:,t+100)=180.*l0/pi;
        ra(:,t+100)=180.*r0/pi;

        %position calculation ---------------------------------------------
        rf(:,4)=lleg(:,7);
        [lleg, lfoot, xcom] = calc_pos(l0,rf,1);
        t0(:,4)=xcom;
        [rleg, rfoot] = calc_rpos(r0,t0);

        lh=cal_lh(lha,t0);
        Rh=cal_rh(rha,t0);
        waist(:,1) = xcom;
        waist(3,1) = waist(3,1) + 25;
        waist(:,2) = xcom;
        waist(:,3) = rleg(:,1);
        waist(:,4) = lleg(:,1);
        zmp = waist(:,2);
        zmp(3)=0;

        
        waist_val(count,1)=waist(3,1);
        count = count + 1;
        
        [z__1, z_1] = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1);
        rlegz(t+100,:) = rleg(3,:);
        llegz(t+100,:) = lleg(3,:);
        pause(wa);
    end
    pause(5);

