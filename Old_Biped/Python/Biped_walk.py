import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time 
from matplotlib import animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
#ax.axis('on')
#ax = fig.add_subplot(111,projection='3d')
#ax.set_xlim(0, 40)
ax.set_ylim(0, 40)
ax.set_zlim(0, 40)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
def Rz(theta):
    r = np.matrix([[np.cos(theta), np.sin(theta), 0, 0], [-np.sin(theta), np.cos(theta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    return r

def cal_lh(lha, t0):
    lh0 = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,25],[0,0,0,1]])*Rz(lha[0])
    lh11 =  np.matrix([[0,0,1,0],[-1,0,0,6],[0,-1,0,0],[0,0,0,1]])*Rz(lha[1])
    lh1 =  np.matrix([[0,0,-1,0],[1,0,0,0],[0,-1,0,0],[0,0,0,1]])*Rz(lha[2])
    lh2 =  np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(lha[3])

    lh[:,0] = t0*np.matrix([[0],[0],[0],[1]])
    lh[:,1] = t0*np.matrix([[0],[0],[25],[1]])
    lh[:,2] = t0*lh0*np.matrix([[0],[0],[5],[1]])
    lh[:,3] = lh[:,1]
    lh[:,4] = t0*lh0*np.matrix([[0],[6],[0],[1]])
    lh[:,5] = t0*lh0*lh11*np.matrix([[0],[0],[0],[1]])
    lh[:,6] = t0*lh0*lh11*lh1*np.matrix([[15],[0],[0],[1]])
    lh[:,7] = t0*lh0*lh11*lh1*lh2*np.matrix([[15],[0],[0],[1]])

    return np.matrix(lh)

def cal_rh(rha,t0):
    rh0 = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,25],[0,0,0,1]])*Rz(rha[0])
    rh11 = np.matrix([[0,0,1,0],[-1,0,0,-6],[0,-1,0,0],[0,0,0,1]])*Rz(rha[1])
    rh1 = np.matrix([[0,0,-1,0],[1,0,0,0],[0,-1,0,0],[0,0,0,1]])*Rz(rha[2])
    rh2 = np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(rha[3])

    Rh[:,0] = t0*np.matrix([[0],[0],[0],[1]])
    Rh[:,1] = t0*np.matrix([[0],[0],[25],[1]])
    Rh[:,2] = t0*rh0*np.matrix([[0],[0],[5],[1]])
    Rh[:,3] = Rh[:,1]
    Rh[:,4] = t0*rh0*np.matrix([[0],[-6],[0],[1]])
    Rh[:,5] = t0*rh0*rh11*np.matrix([[0],[0],[0],[1]])
    Rh[:,6] = t0*rh0*rh11*rh1*np.matrix([[15],[0],[0],[1]])
    Rh[:,7] = t0*rh0*rh11*rh1*rh2*np.matrix([[15],[0],[0],[1]])

    return np.matrix(Rh)

#This function takes 5 inputs (xa, ha, x, h, l) and gives 3 outputs (q3, q4, q5)
def calc_angles(xa, ha, x, h, l):
    x = -x
    q3 = -2.0*np.arctan((4*(l*x - l*xa))/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0) + ((-(h**2.0 - 2.0*h*ha + ha**2.0 + x**2.0 - 2.0*x*xa + xa**2.0)*(h**2.0 - 2.0*h*ha + ha**2.0 - 4*l**2.0 + x**2.0 - 2.0*x*xa + xa**2.0))**(1/2.0) - 2.0*l*x + 2.0*l*xa)/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0))
    q3 = -q3
    q4 = -2.0*np.arctan((4.0*(l*x - l*xa))/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0) + ((-(h**2.0 - 2.0*h*ha + ha**2.0 + x**2.0 - 2.0*x*xa + xa**2.0)*(h**2.0 - 2.0*h*ha + ha**2.0 - 4*l**2.0 + x**2.0 - 2.0*x*xa + xa**2.0))**(1/2.0) - 2.0*l*x + 2.0*l*xa)/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0)) -2.0*np.arctan(((-(h**2.0 - 2.0*h*ha + ha**2.0 + x**2.0 - 2.0*x*xa + xa**2.0)*(h**2.0 - 2.0*h*ha + ha**2.0 - 4*l**2.0 + x**2.0 - 2.0*x*xa + xa**2.0))**(1/2.0) - 2.0*l*x + 2.0*l*xa)/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0))
    q4 = -q4
    q5 = -2.0*np.arctan(((-(h**2.0 - 2.0*h*ha + ha**2.0 + x**2.0 - 2.0*x*xa + xa**2.0)*(h**2.0 - 2.0*h*ha + ha**2.0 - 4*l**2.0 + x**2.0 - 2.0*x*xa + xa**2.0))**(1/2.0) - 2.0*l*x + 2.0*l*xa)/(h**2.0 - 2.0*h*ha + 2.0*l*h + ha**2.0 - 2.0*l*ha + x**2.0 - 2.0*x*xa + xa**2.0))

    return (q3,q4,q5)

def calc_lpos(l0,t0):
    #------------------------ %left leg ----------------------------------------
    lt11 = np.matrix([[0,-1,0,0],[-1,0,0,5],[0,0,-1,0],[0,0,0,1]])*Rz(l0[0])
    lt1 = np.matrix([[0,-1,0,0],[0,0,-1,0],[1,0,0,0],[0,0,0,1]])*Rz(l0[1])
    lt2 = np.matrix([[1,0,0,0],[0,0,1,0],[0,-1,0,0],[0,0,0,1]])*Rz(l0[2])
    lt3 = np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(l0[3])
    lt4 = np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(l0[4])
    lt5 = np.matrix([[1,0,0,0],[0,0,-1,0],[0,1,0,0],[0,0,0,1]])*Rz(l0[5])
    lt6 = np.matrix([[0,1,0,0],[0,0,1,0],[1,0,0,5],[0,0,0,1]])*Rz(l0[6])
    #------------------------- %left leg ---------------------------------------
    lleg[:,0] = t0*lt11*np.matrix([[0],[0],[0],[1]])
    lleg[:,1] = t0*lt11*lt1*np.matrix([[0],[0],[0],[1]])
    lleg[:,2] = t0*lt11*lt1*lt2*np.matrix([[15],[0],[0],[1]])
    lleg[:,3] = t0*lt11*lt1*lt2*lt3*np.matrix([[15],[0],[0],[1]])
    lleg[:,5] = t0*lt11*lt1*lt2*lt3*lt4*lt5*np.matrix([[0],[0],[5],[1]])
    lleg[:,6] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[5],[0],[0],[1]])
    lleg[:,7] = t0*lt11*lt1*lt2*lt3*lt4*lt5*np.matrix([[0],[0],[-5],[1]])
    #------------------------- %left foot --------------------------------------
    lfoot[:,0] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[0],[0],[-4.5],[1]])
    lfoot[:,1] = t0*lt11*lt1*lt2*lt3*lt4*lt5*np.matrix([[0],[-4.5],[-5],[1]])
    lfoot[:,2] = t0*lt11*lt1*lt2*lt3*lt4*lt5*np.matrix([[0],[4.5],[-5],[1]])
    lfoot[:,3] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[0],[0],[4.5],[1]])
    lfoot[:,4] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[5],[0],[4.5],[1]])
    lfoot[:,5] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[5],[0],[-4.5],[1]])
    lfoot[:,6] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[0],[0],[-4.5],[1]])
    lfoot[:,7] = t0*lt11*lt1*lt2*lt3*lt4*lt5*lt6*np.matrix([[0],[0],[4.5],[1]])

    return (lleg, lfoot)

def calc_rpos(r0,t0):
    #right leg --------------------------------------------
    rt11 = np.matrix([[0,-1,0,0],[-1,0,0,-5],[0,0,-1,0],[0,0,0,1]])*Rz(r0[0])
    rt1 = np.matrix([[0,-1,0,0],[0,0,-1,0],[1,0,0,0],[0,0,0,1]])*Rz(r0[1])
    rt2 = np.matrix([[1,0,0,0],[0,0,1,0],[0,-1,0,0],[0,0,0,1]])*Rz(r0[2])
    rt3 = np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(r0[3])
    rt4 = np.matrix([[1,0,0,15],[0,1,0,0],[0,0,1,0],[0,0,0,1]])*Rz(r0[4])
    rt5 = np.matrix([[1,0,0,0],[0,0,-1,0],[0,1,0,0],[0,0,0,1]])*Rz(r0[5])
    rt6 = np.matrix([[0,1,0,0],[0,0,1,0],[1,0,0,5],[0,0,0,1]])*Rz(r0[6])
    #right leg ---------------------------------------------
    rleg[:,0] = t0*rt11*np.matrix([[0],[0],[0],[1]])
    rleg[:,1] = t0*rt11*rt1*np.matrix([[0],[0],[0],[1]])
    rleg[:,2] = t0*rt11*rt1*rt2*np.matrix([[15],[0],[0],[1]])
    rleg[:,3] = t0*rt11*rt1*rt2*rt3*np.matrix([[15],[0],[0],[1]])
    rleg[:,4] = t0*rt11*rt1*rt2*rt3*rt4*np.matrix([[0],[0],[0],[1]])
    rleg[:,5] = t0*rt11*rt1*rt2*rt3*rt4*rt5*np.matrix([[0],[0],[5],[1]])
    rleg[:,6] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[5],[0],[0],[1]])
    rleg[:,7] = t0*rt11*rt1*rt2*rt3*rt4*rt5*np.matrix([[0],[0],[-5],[1]])

    #right foot-------------------------------------------------
    rfoot[:,0] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[0],[0],[-4.5],[1]])
    rfoot[:,1] = t0*rt11*rt1*rt2*rt3*rt4*rt5*np.matrix([[0],[-4.5],[-5],[1]])
    rfoot[:,2] = t0*rt11*rt1*rt2*rt3*rt4*rt5*np.matrix([[0],[4.5],[-5],[1]])
    rfoot[:,3] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[0],[0],[4.5],[1]])
    rfoot[:,4] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[5],[0],[4.5],[1]])
    rfoot[:,5] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[5],[0],[-4.5],[1]])
    rfoot[:,6] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[0],[0],[-4.5],[1]])
    rfoot[:,7] = t0*rt11*rt1*rt2*rt3*rt4*rt5*rt6*np.matrix([[0],[0],[4.5],[1]])

    return (rleg,rfoot)

def calc_pos(r0,t6,wl):
    t5 = t6*((np.linalg.inv(np.matrix([[0,1,0,0],[0,0,1,0],[1,0,0,5],[0,0,0,1]])))*Rz(-r0[5]));
    t4 = t5*((np.linalg.inv(np.matrix([[1 ,0 ,0 ,0],[0 ,0 ,-1 ,0],[0 ,1 ,0 ,0],[0,0,0,1]])))*Rz(-r0[4]));
    t3 = t4*((np.linalg.inv(np.matrix([[1 ,0 ,0 ,15],[0 ,1 ,0 ,0],[0 ,0 ,1 ,0],[0 ,0 ,0 ,1]])))*Rz(-r0[3]));
    t2 = t3*((np.linalg.inv(np.matrix([[1 ,0 ,0 ,15],[0 ,1 ,0 ,0],[0 ,0 ,1 ,0],[0 ,0 ,0 ,1]])))*Rz(-r0[2]));
    t1 = t2*((np.linalg.inv(np.matrix([[1 ,0 ,0 ,0],[0 ,0 ,1 ,0],[0 ,-1 ,0 ,0],[0 ,0 ,0 ,1]])))*Rz(-r0[1]));
    t0 = t1*((np.linalg.inv(np.matrix([[0 ,-1 ,0 ,0],[0 ,0 ,-1 ,0],[1 ,0 ,0 ,0],[0 ,0 ,0 ,1]])))*Rz(-r0[0]));
    leg=np.matrix(np.zeros((4,8)))
    #right leg ---------------------------------------------
    leg[:,0] = t1*np.matrix([[0],[0],[0],[1]]);
    leg[:,1] = t2*np.matrix([[0],[0],[0],[1]]);
    leg[:,2] = t3*np.matrix([[0],[0],[0],[1]]);
    leg[:,3] = t4*np.matrix([[0],[0],[0],[1]]);
    leg[:,4] = t5*np.matrix([[0],[0],[0],[1]]);
    leg[:,6] = t5*np.matrix([[0],[0],[5],[1]]);
    leg[:,7] = t6*np.matrix([[5],[0],[0],[1]]);
    leg[:,5] = t5*np.matrix([[0],[0],[-5],[1]]);

    com = (t0)*([[np.array(wl)*5],[0],[0],[1]]);
    foot=np.matrix(np.zeros((4,8)))
    #right foot-------------------------------------------------
    foot[:,0] = t6*np.matrix([[0],[0],[-4.5],[1]]);
    foot[:,1] = t6*np.matrix([[-10],[0],[-4.5],[1]]);
    foot[:,2] = t6*np.matrix([[-10],[0],[4.5],[1]]);
    foot[:,3] = t6*np.matrix([[0],[0],[4.5],[1]]);
    foot[:,4] = t6*np.matrix([[5],[0],[4.5],[1]]);
    foot[:,5] = t6*np.matrix([[5],[0],[-4.5],[1]]);
    foot[:,6] = t6*np.matrix([[0],[0],[-4.5],[1]]);
    foot[:,7] = t6*np.matrix([[0],[0],[4.5],[1]]);   
    
    return leg,np.matrix(foot),com

def plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1):  ### SOMETHING'S WRONG HERE ###
    ax.clear()
    #plt.hold(False)
    #print np.array(rleg[2,:]).T
    #ax.plot3D(np.array(rleg[0,:]).T,np.array(rleg[1,:]).T,np.array(rleg[2,:]).T,'-o')
    #ax.plot(np.array(lleg[0,:]).T,np.array(lleg[1,:]).T,np.array(lleg[2,:]).T,'k-')
    #ax.scatter3D(np.array(zmp[0]).T,np.array(zmp[1]).T,np.array(zmp[2]).T,c=np.array(zmp[2]).T, cmap='Greens')
    #ax.plot(np.array(rfoot[0,:]).T,np.array(rfoot[1,:]).T,np.array(rfoot[2,:]).T,'blue')
    #ax.plot(np.array(lfoot[0,:]).T,np.array(lfoot[1,:]).T,np.array(lfoot[2,:]).T,'red')        
    #ax.plot3D(waist[0,:],waist[1,:],waist[2,:],'-o')
    #ax.plot3D(zmp[0],zmp[1],zmp[2],'-o')
    #ax.scatter(np.array(zmp[0]).T,np.array(zmp[1]).T,np.array(zmp[2]).T)
    ax.scatter(np.array(lleg[0,:]).T,np.array(lleg[1,:]).T,np.array(lleg[2,:]).T)
    ax.scatter(np.array(rleg[0,:]).T,np.array(rleg[1,:]).T,np.array(rleg[2,:]).T)
    ax.scatter(np.array(zmp[0]).T,np.array(zmp[1]).T,np.array(zmp[2]).T)
    ax.scatter(np.array(rfoot[0,:]).T,np.array(rfoot[1,:]).T,np.array(rfoot[2,:]).T)
    ax.scatter(np.array(lfoot[0,:]).T,np.array(lfoot[1,:]).T,np.array(lfoot[2,:]).T)
    ax.scatter(np.array(waist[0,:]).T,np.array(waist[1,:]).T,np.array(waist[2,:]).T)
    ax.plot_wireframe(np.array(lleg[0,:]).T,np.array(lleg[1,:]).T,np.array(lleg[2,:]).T)
    ax.plot_wireframe(np.array(rleg[0,:]).T,np.array(rleg[1,:]).T,np.array(rleg[2,:]).T)
    ax.plot_wireframe(np.array(zmp[0]).T,np.array(zmp[1]).T,np.array(zmp[2]).T)
    ax.plot_wireframe(np.array(rfoot[0,:]).T,np.array(rfoot[1,:]).T,np.array(rfoot[2,:]).T)
    ax.plot_wireframe(np.array(lfoot[0,:]).T,np.array(lfoot[1,:]).T,np.array(lfoot[2,:]).T)
    ax.plot_wireframe(np.array(waist[0,:]).T,np.array(waist[1,:]).T,np.array(waist[2,:]).T)
    #ax.plot3D(rfoot[0,:],rfoot[1,:],rfoot[2,:],'-o')
   # ax.plot3D(lfoot[0,:],lfoot[1,:],lfoot[2,:],'-o')
    #ax.set_xlim(-40, 40)
    ax.set_ylim(-40, 40)
    ax.set_zlim(0, 40)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    #print(rleg)
    #print(lleg)
    #print(waist)
    #print(zmp)
    #print(rfoot)
    #print(lfoot)
    #plt.show()
           
    #plot3(rleg(1,:),rleg(2,:),rleg(3,:),'k -o',lleg(1,:),lleg(2,:),lleg(3,:),'k -o',waist(1,:),waist(2,:),waist(3,:),'k -o',zmp(1),zmp(2),zmp(3),'b -o',rfoot(1,:),rfoot(2,:),rfoot(3,:),'k -o',lfoot(1,:),lfoot(2,:),lfoot(3,:),'k -o','LineWidth',1.5);

    #view([0 0]);
    z=30*np.cos(np.pi*15/180);
    d=z/(980*0.000001);
    a=-d;
    b=2*d+1;
    c=-d;

    p = a*z__1 + b*z_1 + c*zmp;
    #plot(z_1(1),z_1(2),'b--o');
    #plot(zmp(1),zmp(2),'b--o');
    #if (p(2)<10) && (p(2)>-10)
     #   plot(p(1),p(2),'r--o');
    #end
    
    #grid on;
    #box on;
    #axis equal;
    plt.pause(0.005)
    #axis square;
    #axis([-15 150 -15 15 -15 50]);
    #hold off;
    z2 = z_1;
    z1 = zmp;

    return (z2,z1)


##### Variables Declarartion of appropriate size ##############

Rh=np.matrix(np.zeros((4,8)))
lh=np.matrix(np.zeros((4,8)))
rha =np.matrix(np.zeros((4,1)))
lha =np.matrix(np.zeros((4,1)))
waist=np.matrix(np.zeros((4,4)))
rleg=np.matrix(np.zeros((4,8)))
lleg=np.matrix(np.zeros((4,8)))
rlegz=np.matrix(np.zeros((200,8)))
llegz=np.matrix(np.zeros((200,8)))
l0=np.matrix(np.zeros((7,1)))
r0=np.matrix(np.zeros((7,1)))
rfoot=np.matrix(np.zeros((4,8)))
lfoot=np.matrix(np.zeros((4,8)))
xcom=np.matrix(np.zeros((4,1)))
zmp =np.matrix(np.zeros((4,1)))
z=30*np.cos(np.pi*18/180)
h=z
ll = 15.0
y=0.0
t=0
t0 = np.matrix([[1.0, 0, 0, 0],[0, 1.0, 0, 0],[0, 0, 1.0, 0],[0, 0, 0, 1.0]])
rf = np.matrix([[1.0, 0, 0, 0],[0, 0, 1.0, 0],[0, -1.0, 0, 0],[0, 0, 0, 1.0]])
wa1=0.001
wa2=0.001
wa = 0.01
sl=16.0 #step length
z__1=np.matrix(np.zeros((4,1)))
z_1=np.matrix(np.zeros((4,1)))
la=np.matrix(np.zeros((7,200)))
ra=np.matrix(np.zeros((7,200)))




    ###############################################################
    ######    right leg as swing leg                          #####
    ###############################################################

for t in range(1,100):
    #t=0

    print t
    x=t/10.0

    if (t>-1) and (t<41): # SSP1  (or DSP1 ?????)
        y=-(3*t)/40.0
        xa = -5.0
        ha = 0.0
    elif (t>40) and (t<81): # DSP1  (or SSP1 ?????)
        y = -3.0
        xa = (20/4)*(x-4.0) - 5
        ha = 0.00237671*(xa)**3 -0.0582294*(xa)**2 + 0.16637*(xa) + 2.58467
    else: # SSP2  (or DSP2 ?????)
        y = 3*(t-80)/20.0 - 3
        xa = 15.0
        ha = 0



    l0[5]=np.arctan((-y)/h)
    l0[1]=-l0[5]
    (l0[2], l0[3], l0[4]) = calc_angles(-xa,ha,x,h,ll)
    l0[4]=-l0[4]
    l0[3]=-l0[3]

    r0[5] = l0[5]
    r0[1] = l0[1]
    r0[2], r0[3], r0[4] = calc_angles(-5,0,x,h,ll)
    r0[3]=-r0[3]
    r0[4]=-r0[4]

    la[:,t-1]=180.*l0/np.pi
    ra[:,t-1]=180.*r0/np.pi

    #position calculation ---------------------------------------------
    rf[:,3]=rleg[:,6]
    rleg, rfoot, xcom = calc_pos(r0,rf,-1)
    t0[:,3]=xcom
    lleg, lfoot = calc_lpos(l0,t0)
    lh=cal_lh(lha,t0)
    Rh=cal_rh(rha,t0)
    waist[:,0] = xcom
    waist[2,0] = waist[2,0] + 25
    waist[:,1] = xcom
    waist[:,2] = rleg[:,0]
    waist[:,3] = lleg[:,0]
    zmp = np.array(waist[:,1])
    zmp=np.matrix(zmp)
    zmp[2]=0
    z__1, z_1 = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1)
    rlegz[t-1,:] = rleg[2,:]
    llegz[t-1,:] = lleg[2,:]
    #time.sleep(wa) #for this just add t.sleep(wa) and then comment it.


###############################################################
######    left leg as swing leg                          #####
###############################################################

for t in range(1,100):
    #t=0;
    x=(t)/10.00;
    if (t>-1) and (t<41): # SSP1  (or DSP1 ?????)
        y=(3*t)/40.0
        xa = -5.0
        ha = 0.0
    elif (t>40) and (t<81): # DSP1  (or SSP1 ?????)
        y = 3.0
        xa = (20/4)*(x-4.0) - 5
        ha = 0.00237671*(xa)**3 -0.0582294*(xa)**2 + 0.16637*(xa) + 2.58467
    else: # SSP2  (or DSP2 ?????)
        y = -3*(t-80)/20.0 - 3
        xa = 15.0
        ha = 0
    l0[5]=np.arctan((-y)/h);
    l0[1]=-l0[5];
    l0[2], l0[3], l0[4] = calc_angles(-xa,ha,x,h,ll);
    l0[4]=-l0[4];
    l0[3]=-l0[3];
    r0[5] = l0[5];
    r0[1] = l0[1];
    r0[2], r0[3], r0[4] = calc_angles(-5,0,x,h,ll);
    r0[3]=-r0[3];
    r0[4]=-r0[4];

    r0,l0=l0,r0  #swap

    la[:,t+99]=180.*l0/np.pi;
    ra[:,t+99]=180.*r0/np.pi;

    #position calculation ---------------------------------------------
    rf[:,3]=lleg[:,6];
    lleg, lfoot, xcom = calc_pos(l0,rf,1);
    t0[:,3]=xcom;
    rleg, rfoot = calc_rpos(r0,t0);

    lh=cal_lh(lha,t0);
    Rh=cal_rh(rha,t0);
    waist[:,0] = xcom;
    waist[2,0] = waist[2,0] + 25;
    waist[:,1] = xcom;
    waist[:,2] = rleg[:,0];
    waist[:,3] = lleg[:,0];
    zmp = np.array(waist[:,1])
    zmp=np.matrix(zmp)
    zmp[2]=0;

    z__1, z_1 = plot_final(Rh,lh,rleg,lleg,waist,zmp,rfoot,lfoot,z__1,z_1);
    
    rlegz[t+99,:] = rleg[2,:]
    llegz[t+99,:] = lleg[2,:]
    #time.sleep(wa)