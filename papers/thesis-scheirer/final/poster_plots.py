from __future__ import division
import numpy as np
import sys
import matplotlib as mpl
import SW
import pylab as plt



from matplotlib.path import Path
from matplotlib.patches import PathPatch

data=np.loadtxt('figs/snaft2.out')

x = plt.linspace(1e-20,.17,40000)                         # x-axis grid space (used for when plotting number density on x-axis)

xmin=((4*np.pi*(SW.R)**3)/3)*1e-20
xmax=((4*np.pi*(SW.R)**3)/3)*.17
#xmin=(3/(4.0*np.pi))*(1e-20)*(1/(SW.R)**3)
#xmax=(3/(4.0*np.pi))*(0.2)*(1/(SW.R)**3)
xff = plt.linspace(xmin,xmax,40000)








nL = [data[i][1] for i in range(0,len(data))]           # list of the left number density solutions obtained from fsolve
nR = [data[i][2] for i in range(0,len(data))]           # list of the right number density solutions obtained from fsolve
Tlist = [data[i][0] for i in range(0,len(data))]        # list of the corresponding temperatures for which the above number densitites were found

ffL = [i*((4*np.pi*(SW.R)**3)/3) for i in nL]           # converts left number density to filling fraction
ffR = [i*((4*np.pi*(SW.R)**3)/3) for i in nR]           # converts right number density to filling fraction



def liq_vap_Tvsff():
  plt.figure()
  plt.plot(ffL,Tlist,color='#f36118',linewidth=2)
  plt.plot(ffR,Tlist,color='c',linewidth=2)
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  plt.xlim(-.05,max(ffR)+.05)
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
  #plt.savefig('figs/liqVapCo_Tvsff.pdf')


  plt.figure()
  plt.plot(ffL,Tlist,color='#f36118',linewidth=2)
  plt.plot(ffR,Tlist,color='c',linewidth=2)
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  plt.xlim(-.05,max(ffR)+.05)
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
  fig=plt.figure(figsize=(1,5))
  ax=fig.add_subplot(111)
  ax.axis([0,1,-50,200])
  cmap = mpl.cm.jet
  norm = mpl.colors.Normalize(vmin=-40, vmax=180)
  cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                orientation='vertical',                                
                                norm=norm,
                                ticks=[-40,180]
                                )

  plt.subplots_adjust(left=0.4, right=0.8)
  
  #plt.savefig('figs/liqVapCo_Tvsff2.pdf')
 

  plt.figure()
  xx=np.arange(0,10,0.01)
  yy=xx*np.exp(-xx)

  path = Path(np.array([xx,yy]).transpose())
  patch = PathPatch(path, facecolor='none')
  plt.gca().add_patch(patch)

  im = plt.imshow(xx.reshape(yy.size,1),  cmap=plt.cm.Reds,interpolation="bicubic",
                origin='lower',extent=[0,10,-0.0,0.40],aspect="auto", clip_path=patch, clip_on=True)
  #im.set_clip_path(patch)
  #plt.savefig("out.png")

  plt.show()



  

def liq_vap_co_pvsT():
        pL=[]
        pR=[]
        pdiff=[]
        for i in range(0,len(nL)):
                pL.append(SW.findP(Tlist[i],nL[i]))
                pR.append(SW.findP(Tlist[i],nR[i]))
                pdiff.append(pL[i]-pR[i])
        plt.figure()
        plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
        plt.ylabel('pressure')
        plt.xlabel('temperature')
        plt.plot(Tlist,pL,color='#f36118',linewidth=8)
        plt.plot(Tlist,pR,'c',linewidth=3)
        plt.savefig('figs/liqVapCo_pvsT.pdf')
        #plt.show()
        
        #plt.figure()
        #plt.title('Pressure check')
        #plt.ylabel('P')
        #plt.xlabel('T')
        #plt.plot(Tlist,pdiff,'r')
        #plt.savefig('figs/liqVapCo_pvsT_diff.pdf')


def cotangent_t0():
        plt.figure()
        plt.title('Helmholtz energy per volume  VS  ff  @ T=%0.4f'%Tlist[100])
        plt.ylabel('Helmholtz free energy per volume')
        plt.xlabel('filling fraction')
        plt.plot(xff,SW.ftot(Tlist[100],x),color='#f36118',linewidth=3)
        plt.plot(xff, SW.numH_dftot_dn(Tlist[100],nR[100])*(x-nR[100])+SW.ftot(Tlist[100],nR[100]),color='#00c0c0',linewidth=2)
        #plt.plot(nL[100],SW.ftot(Tlist[100],nL[100]),'ko')
        #plt.plot(nR[100],SW.ftot(Tlist[100],nR[100]),'ko')
        plt.plot(nL[100]*((4*np.pi*(SW.R)**3)/3),SW.ftot(Tlist[100],nL[100]),'ko')
        plt.plot(nR[100]*((4*np.pi*(SW.R)**3)/3),SW.ftot(Tlist[100],nR[100]),'ko')
        plt.savefig('figs/hfe_cotangent.pdf')
        #plt.show()


def gfe():
  x2=plt.linspace(1e-20,.12,40000)
  mu=SW.numH_dftot_dn(Tlist[100],nR[100])
  plt.figure()
  plt.title('Grand free energy per volume vs n@ T=%0.4f'%Tlist[100])
  plt.ylabel('Grand free energy per volume')
  plt.xlabel('number density (n)')  
  plt.plot(x2,SW.ftot(Tlist[100],x2)-mu*x2,color='#f36118')
  plt.plot(x2,SW.numH_dftot_dn(Tlist[100],nR[100])*(x2-nR[100])+SW.ftot(Tlist[100],nR[100])-mu*nR[100],'c')
  plt.plot(nL[100],SW.ftot(Tlist[100],nL[100])-mu*nR[100],'ko')
  plt.plot(nR[100],SW.ftot(Tlist[100],nR[100])-mu*nR[100],'ko')
  plt.show()


def gfe2():
  x2=plt.linspace(1e-20,.12,40000)
  mu=SW.numH_dftot_dn(Tlist[100],nR[100])*(x2-nR[100])+SW.ftot(Tlist[100],nR[100])
  plt.figure()
  plt.title('Grand free energy per volume vs n@ T=%0.4f'%Tlist[100])
  plt.ylabel('Grand free energy per volume')
  plt.xlabel('number density (n)')  
  plt.plot(x2,SW.ftot(Tlist[100],x2)-mu,color='#f36118')
  plt.plot(x2,x2-x2,'c')
  plt.plot(nL[100],0,'ko')
  plt.plot(nR[100],0,'ko')
  
  #plt.savefig('figs/gfe2.pdf')
  plt.show()
  
def gfe3():
  x2=plt.linspace(1e-20,.2,40000)
  mu=SW.numH_dftot_dn(Tlist[100],x2)
  plt.figure()
  plt.title('Grand free energy per volume vs n@ T=%0.4f'%Tlist[100])
  plt.ylabel('Grand free energy per volume')
  plt.xlabel('number density (n)')  
  plt.plot(x2,SW.ftot(Tlist[100],x2)-mu*nR[100],color='#f36118')
  plt.plot(x2,x2-x2,'c')
  plt.plot(nL[100],0,'ko')
  plt.plot(nR[100],0,'ko')
  plt.show()


def gfe4():
  x2=plt.linspace(1e-20,.13,90000)
  xmin2=((4*np.pi*(SW.R)**3)/3)*1e-20
  xmax2=((4*np.pi*(SW.R)**3)/3)*.13
  xff2 = plt.linspace(xmin2,xmax2,90000)
  thigh=100
  plt.figure()
  plt.title('Grand free energy per volume vs ff @ T=%0.4f'%Tlist[thigh])
  plt.ylabel('Grand free energy per volume')
  plt.xlabel('filling fraction')  
  plt.plot(xff2,SW.phi(Tlist[thigh],x2,nR[thigh]),color='#f36118',linewidth=3)
  #plt.axvline(nL[thigh])
  #plt.axvline(nR[thigh])
  #plt.axhline(SW.phi(Tlist[thigh],nR[thigh]))
  #plt.plot(x2,x2-x2,'c')
  plt.plot(nL[thigh]*((4*np.pi*(SW.R)**3)/3),SW.phi(Tlist[thigh],nL[thigh],nR[thigh]),'ko')
  plt.plot(nR[thigh]*((4*np.pi*(SW.R)**3)/3),SW.phi(Tlist[thigh],nR[thigh],nR[thigh]),'ko')
  plt.axhline(SW.phi(Tlist[thigh],nR[thigh],nR[thigh]),color='c',linewidth=2)
  print(Tlist[100])
  print(nL[100],nR[100])
  plt.savefig('figs/gfe_cotangent.pdf')

  plt.figure()
  plt.plot(xff2,SW.phi(Tlist[thigh],x2,nR[thigh]),color='#f36118',linewidth=3)
  plt.plot(nL[thigh]*((4*np.pi*(SW.R)**3)/3),SW.phi(Tlist[thigh],nL[thigh],nR[thigh]),'ko')
  plt.plot(nR[thigh]*((4*np.pi*(SW.R)**3)/3),SW.phi(Tlist[thigh],nR[thigh],nR[thigh]),'ko')
  plt.axhline(SW.phi(Tlist[thigh],nR[thigh],nR[thigh]),color='c',linewidth=2)
  plt.xlim(0,0.0003)
  plt.ylim(-.000014,0.000006)
  print(Tlist[100])
  print(nL[100],nR[100])
  plt.savefig('figs/gfe_insert_cotangent.pdf')
  
  #plt.show()


def ideal():
  x2=plt.linspace(1e-20,.13,40000)
  
  tt=plt.linspace(1e-20,1.2,4000)
  plt.figure()
  plt.title('ideal gas')
  plt.ylabel('f')
  plt.xlabel('n')  
  #plt.plot(x2,SW.fid(Tlist[100],x2),color='#f36118')
  #plt.plot(tt,SW.fid(tt,nR[100]))
  plt.show()


def gfe5():
  x2=plt.linspace(1e-20,.2,4000)
  lazy=100
  mu=SW.numH_dftot_dn(Tlist[lazy],x2)
  plt.figure()
  plt.title('Grand free energy per volume vs n@ T=%0.4f'%Tlist[lazy])
  plt.ylabel('Grand free energy per volume')
  plt.xlabel('number density (n)')  
  plt.plot(x2,SW.ftot(Tlist[lazy],x2)-mu*x2,color='#f36118')
  plt.plot(x2,x2-x2,'c')
  plt.plot(nL[lazy],0,'ko')
  plt.plot(nR[lazy],0,'ko')
  plt.show()


def liq_vap_Tvsn():
  plt.figure()
  plt.plot(nL,Tlist,color='#f36118',linewidth=2)
  plt.plot(nR,Tlist,color='c',linewidth=2)
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  plt.xlim(-.05,max(nR)+.05)
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
  #plt.savefig('figs/liqVapCo_Tvsff.pdf')


  '''plt.figure()
  plt.plot(nL,Tlist,color='#f36118',linewidth=2)
  plt.plot(nR,Tlist,color='c',linewidth=2)
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  plt.xlim(-.05,max(nR)+.05)
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')'''

  plt.show()


def cotangent_t100():
        x2=plt.linspace(1e-20,.2,4000)
        plt.figure()
        plt.title('Helmholtz energy per volume  VS  ff  @ T=%0.4f'%Tlist[100])
        plt.ylabel('Helmholtz free energy per volume')
        plt.xlabel('filling fraction')
        plt.plot(x2,SW.ftot(Tlist[100],x2),color='#f36118',linewidth=3)
        plt.plot(x2, SW.numH_dftot_dn(Tlist[100],nR[100])*(x2-nR[100])+SW.ftot(Tlist[100],nR[100]),color='#00c0c0',linewidth=2)
        #plt.plot(nL[100],SW.ftot(Tlist[100],nL[100]),'ko')
        #plt.plot(nR[100],SW.ftot(Tlist[100],nR[100]),'ko')
        plt.plot(nL[100],SW.ftot(Tlist[100],nL[100]),'ko')
        plt.plot(nR[100],SW.ftot(Tlist[100],nR[100]),'ko')
        #plt.savefig('figs/hfe_cotangent.pdf')
        plt.show()

#ideal()
#liq_vap_co_pvsT()
#liq_vap_Tvsff()
#cotangent_t0()
#gfe()
#gfe2()
#gfe3()
gfe4()
#liq_vap_Tvsn()
#cotangent_t100()
        
