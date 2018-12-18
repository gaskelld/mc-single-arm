      subroutine xsec_model(ispec,ebeam,pcent,thcent,delta,yptar,xptar,
     > xbj,xsec)

      implicit none

      real*8 ebeam,pcent,thcent
      real*8 delta,yptar,xptar
      real*8 plabx,plaby,plabz,pz,p
      real*8 costheta,sintheta
      real*8 thetalab,raddeg
      real*8 nu,wsq,q2,xbj,mp,qabs
      real*8 A, Z, m_tgt
      real*8 sig_dis, sig_qe,sigtot
      real*8 xsec
      integer ispec

      mp=0.938272
      raddeg = 180.0/3.14159

      p=(pcent/1000.0)*(1.0+delta/100.0)
      pz=p/sqrt(1.0+xptar**2+yptar**2)

      sintheta=sin(thcent)
      costheta=cos(thcent)
      if(ispec.eq.2) then
         plabx=pz*xptar
         plaby=pz*(yptar*costheta+sintheta)
         plabz=pz*(-yptar*sintheta+costheta)
      elseif(ispec.eq.1) then
         plabx=pz*xptar
         plaby=pz*(yptar*costheta-sintheta)
         plabz=pz*(yptar*sintheta+costheta)
      endif

      thetalab=acos(plabz/p)

      qabs=sqrt((plabx**2+plaby**2+(ebeam-plabz)**2))
c      write(6,*) 'bad kitty',plabx,plaby,plabz
      nu=Ebeam-p
      q2=qabs**2-nu**2
      wsq=-q2+mp**2+2.*mp*nu
      xbj=q2/2./mp/nu
c      write(6,*) 'cheesy poofs',thetalab*raddeg,nu,q2,xbj

C Carbon target, 1.5%
      A=12.0
      Z=6.0
      m_tgt=12.0*0.931494
      call sigmodel_calc(ebeam,p,thetalab*raddeg,z,a,m_tgt,sig_dis,
     >              sig_qe,sigtot) !nb/sr/MeV


      xsec=sigtot
      return
      end

