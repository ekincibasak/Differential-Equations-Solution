!//This  fortran code  calculates the flavor evolution of neutrinos in the Supernova.

!nonlinear differential equation by fourth order runge-kutta
!// Input:
!// Energy vector E. It should be given in MeV
!// r: vector of radius values. 
!// Outputs:
!// rho: The density operatorof neutrinos.
!// neutrino flux on the end of the supernovae
!Runge-Kutta driver with adaptive stepsize control. Integrate the array of starting values
!ystart (psiim) from x1(r_begin) to x2(r_final) with accuracy eps , storing intermediate results in the module
!variables in ode path . h1 should be set as a guessed first stepsize, hmin as the minimum
!allowed stepsize (can be zero). On output ystart is replaced by values at the end of the
!integration interval. derivs is the user-supplied subroutine for calculating the right-hand-
!side derivative, while rkqc is the name of the stepper routine to be used.
 Program Tequdif
  implicit none

  integer, parameter :: ned=500, nmd=3
  integer kdn,nin,nvar
  real(kind=kind(8.0D0))  h,t1,t2,Enn,numb,L_g,F2,nbad,nok,R,rcs
  integer end_val,intme(nmd,nmd) 
  complex(kind=kind(8.0D0)) psiim(ned,nmd,nmd),Vnpadm(nmd,nmd)
  complex(kind=kind(8.0D0)) apsiim(ned,nmd,nmd) 
  real(kind=kind(8.0D0)) Tme(nmd,nmd),Eme(nmd,nmd)

  integer :: i,j, nmax
  double precision, dimension(9,500) ::  datei
  double precision, dimension(10000000) :: xxx,yyy
  real, parameter :: pi=3.141592654

  real(kind=kind(8.0D0))  aEe,aEm,aEt,Ee,Em,Et,aTm,aTt,aTe,Tm,Tt,Te
  real*8 nse,nsm,nst,anse,ansm,anst,g
  real*8 erho,mrho,trho,aerho,amrho,atrho,TTN,ATTN
  real*8 ndfe,ndfm,ndft,andfe,andfm,andft 
  real*8 rhoce,rhocm,rhoct,arhoce,arhocm,arhoct

  OPEN (UNIT=7, FILE="IH_5in_disneut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=8, FILE="IH_5anin_disneut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=9, FILE="IH_5in_neut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=10, FILE="IH_5anin_neut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=14, FILE="IH_5fi_disneut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=15, FILE="IH_5anfi_disneut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=16, FILE="IH_5neut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=17, FILE="IH_5anneut.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=18, FILE="IH_5neutfidist.dat", ACTION="write", STATUS="unknown")
  OPEN (UNIT=19, FILE="IH_5anneutfidist.dat", ACTION="write", STATUS="unknown")

  open(1, file='shockdensity2.dat',status='old')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  read  (*,*)  h !stepsize
  read  (*,*) t1,t2 !starting time(t1), ending time(t2)
  read  (*,*) nvar !dimension of matrix
  read  (*,*) end_val,nin !to adjust energy interval
  read  (*,*) numb 
!!!Self interaction parameters
  read  (*,*) Te,Tm,Tt  !Te,Tm,Tt tempture
  read  (*,*) aTe,aTm,aTt  !aTe,aTm,aTt
  read  (*,*) R !Radius
  read  (*,*) L_g,g  !luminosity (MeV/km), fermiintegral of order


!Tme matris include neutrino, antineutrino tempture values.
  Tme(1,1)=Te
  Tme(1,2)=Tm
  Tme(1,3)=Tt
  Tme(2,1)=aTe
  Tme(2,2)=aTm
  Tme(2,3)=aTt
  Ee=(3.1514*Tme(1,1)); ;!electron neutrino temperature, Mev
  Em=(3.1514*Tme(1,2)); !muon neutrino temperature, Mev
  Et=(3.1514*Tme(1,3)); !tau neutrino temperature, Mev
  aEe=(3.1514*Tme(2,1)) !electron neutrino temperature, Mev
  aEm=(3.1514*Tme(2,2)); !muon neutrino temperature, Mev
  aEt=(3.1514*Tme(2,3)); !tau neutrino temperature, Mev


!Eme matris include neutrino, antineutrino energy values and R,L_g and fermi couple constant.
  Eme(1,1)=Ee
  Eme(1,2)=Em
  Eme(1,3)=Et

  Eme(2,1)=aEe
  Eme(2,2)=aEm
  Eme(2,3)=aEt

  Eme(3,1)=R
  Eme(3,2)=L_g
  Eme(3,3)=g !alpha+1 =pincf parameter. !alpha=3 !Keil_03
!intme matris include the inial and final values of neutrino energy.
  intme(1,1)=end_val
  intme(1,2)=nin
  !intme(1,3)=F2
!rcs 
  rcs=1/(4*pi*(Eme(3,1)**2)) !

  Tme(3,1)=rcs

! read the density(gr/cm^3) and radius (km)
  do i=1,500
     read(1,*,end=11) (datei(i,j),j=1,9)
     !print*, i
     xxx(i)= datei(i,1) !km
     yyy(i)= datei(i,6) !gr/cm^3
     nmax=i
  end do
11 close(1)	

  do 10 i=1,intme(1,2)
     kdn=0+i

!define the neutrino energy range!0.5-100.5
     Enn = (intme(1,1)*(i+1)*0.5)/intme(1,2)
!the normalized neutrino spectra
     nse= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(1,1)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(1,1));   !!!!loook at!!!!!!!
     nsm= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(1,2)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(1,2));
     nst= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(1,3)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(1,3));
     anse= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(2,1)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(2,1));
     ansm= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(2,2)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(2,2));
     anst= (((Eme(3,3)**Eme(3,3))/6)*((Enn**(Eme(3,3)-1))/(Eme(2,3)**Eme(3,3)) ))&
          &* exp(-(Eme(3,3)*Enn)/Eme(2,3));


!Em(3,2) is the Luminosity
! Luminosity [MeV/km] => (1/6)*1e+53 erg/s ;(yoshida_04)
! 1 [erg/s] = 2.084 [MeV/km]
! T_supernova => For Shockwave
!the number density of specific neutrinos of pencil direction. (but(1/piR^2) and Dr part is missing.
!Unit is 1/MeVkm
     rhoce=(Eme(3,2)*(1/(Eme(1,1)))*nse)
     rhocm=(Eme(3,2)*(1/(Eme(1,2)))*nsm)
     rhoct=(Eme(3,2)*(1/(Eme(1,3)))*nst)



     arhoce=(Eme(3,2)*(1/(Eme(2,1)))*anse)
     arhocm=(Eme(3,2)*(1/(Eme(2,2)))*ansm)
     arhoct=(Eme(3,2)*(1/(Eme(2,3)))*anst)




!It is the trace of the matrix give the total number.
     TTN=rhoce+rhocm+rhoct
     !en*ede+mn*edm+tn*edt
     ATTN=arhoce+arhocm+arhoct

     !aen*aede+amn*aedm+atn*aedt
!ratio of the specic neutrino number density to total number density.
     erho=rhoce/TTN
     mrho=rhocm/TTN
     trho=rhoct/TTN


     aerho=arhoce/ATTN
     amrho=arhocm/ATTN
     atrho=arhoct/ATTN




     psiim(kdn,1,1)=erho
     psiim(kdn,1,2)=0.d0
     psiim(kdn,1,3)=0.d0
     psiim(kdn,2,1)=0.d0
     psiim(kdn,2,2)=mrho; 
     psiim(kdn,2,3)=0.d0
     psiim(kdn,3,1)=0.d0
     psiim(kdn,3,2)=0.d0
     psiim(kdn,3,3)=trho

     apsiim(kdn,1,1)=aerho
     apsiim(kdn,1,2)=0.d0
     apsiim(kdn,1,3)=0.d0
     apsiim(kdn,2,1)=0.d0
     apsiim(kdn,2,2)=amrho; 
     apsiim(kdn,2,3)=0.d0
     apsiim(kdn,3,1)=0.d0
     apsiim(kdn,3,2)=0.d0
     apsiim(kdn,3,3)=atrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!rcs=1/(4*pi*(Eme(3,1)**2))1/km

!Tme(3,1)=rcs
!unit is MeV/km^3 to MeV^2  multriply 1.97*10^-16

     WRITE (UNIT=7, FMT=*) Enn,t1,real(psiim(kdn,1,1))*(TTN*Tme(3,1)),real(psiim(kdn,1,2))*TTN*Tme(3,1),&
          &real(psiim(kdn,1,3))*TTN*Tme(3,1),real(psiim(kdn,2,1))*TTN*Tme(3,1),real(psiim(kdn,2,2))*(TTN*Tme(3,1)),&
          &real(psiim(kdn,2,3))*TTN*Tme(3,1),real(psiim(kdn,3,1))*TTN*Tme(3,1),real(psiim(kdn,3,2))*TTN*Tme(3,1),&
          &real(psiim(kdn,3,3))*(TTN*Tme(3,1))

     WRITE (UNIT=8, FMT=*)Enn,t1,real(apsiim(kdn,1,1))*(ATTN*Tme(3,1)),real(apsiim(kdn,1,2))*(ATTN*Tme(3,1)),&
          &real(apsiim(kdn,1,3))*(ATTN*Tme(3,1)),real(apsiim(kdn,2,1))*(ATTN*Tme(3,1)),real(apsiim(kdn,2,2))*(ATTN*Tme(3,1)),&
          &real(apsiim(kdn,2,3))*(ATTN*Tme(3,1)),real(apsiim(kdn,3,1))*(ATTN*Tme(3,1)),real(apsiim(kdn,3,2))*(ATTN*Tme(3,1)),&
          &real(apsiim(kdn,3,3))*(ATTN*Tme(3,1))



     WRITE (UNIT=9, FMT=*)Enn,t1,real(psiim(kdn,1,1)),real(psiim(kdn,1,2)),&
          &real(psiim(kdn,1,3)),real(psiim(kdn,2,1)),real(psiim(kdn,2,2)),&
          &real(psiim(kdn,2,3)),real(psiim(kdn,3,1)),real(psiim(kdn,3,2)),&
          &real(psiim(kdn,3,3))

     WRITE (UNIT=10, FMT=*)Enn,t1,real(apsiim(kdn,1,1)),real(apsiim(kdn,1,2)),&
          &real(apsiim(kdn,1,3)),real(apsiim(kdn,2,1)),real(apsiim(kdn,2,2)),&
          &real(apsiim(kdn,2,3)),real(apsiim(kdn,3,1)),real(apsiim(kdn,3,2)),&
          &real(apsiim(kdn,3,3))

10   continue 


     call odeint(psiim,apsiim,nvar,numb,xxx,yyy,nmax,t1 ,t2,1.d-6,h,0.000000000001d0,nok,nbad,&
          &Tme,Eme,intme,Vnpadm)

  END do

  Subroutine Derivs(yy,rtt,psim,apsim,numb,dxxf,adxxf,Tme,Eme,intme,Vnpadm)
!Subroutine Derivs incluede both the derivatives at the starting point and the function. 

    real(kind=kind(8.0D0)) Dr,g,numb,const

    !real*8, parameter :: tmax=1000.d0

    complex*8,parameter :: ic=( 0.0, 1 )
    real, parameter :: pi=3.141592654
    double precision :: yy,rtt,yyh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: ned=500, nmd=3
    integer intme(nmd,nmd)



    real(kind=kind(8.0D0)) Tme(nmd,nmd),Eme(nmd,nmd)
    real(kind=kind(8.0D0)) ::Hvv(nmd,nmd)
    complex(kind=kind(8.0D0)) dxxf(ned,nmd,nmd),adxxf(ned,nmd,nmd),Vnp(nmd,nmd),sint(nmd,nmd)      
    integer :: c1(nmd,nmd)
    complex(kind=kind(8.0D0)) psim(ned,nmd,nmd),apsim(ned,nmd,nmd),dxxfn(nmd,nmd),adxxfn(nmd,nmd) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    real(kind=kind(8.0D0)) R,L_g,h,radfc
    real(kind=kind(8.0D0)) TTN,ATTN,eigenvalue,eigenvaluea,element12

    real(kind=kind(8.0D0))omega,delta_m_sqr31,theta13,omegcos,omegsinsqr
    real(kind=kind(8.0D0)) eigenval,theta,tht

    real*8 aEe,aEm,aEt,Ee,Em,Et,aTm,aTt,aTe,Tm,Tt,Te
    real*8 nse,nsm,nst,anse,ansm,anst
    real*8 rhoce,rhocm,rhoct,arhoce,arhocm,arhoct



    real(kind=kind(8.0D0)) rhocek,rhocmk,rhoctk,arhocek,arhocmk,arhoctk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: nin,nkk,nkkdz

    real :: u,hint


    complex(kind=kind(8.0D0)) integral(nmd,nmd),valueadm(nmd,nmd),&
         &kselfin,valuek(nmd,nmd),integraladm(nmd,nmd),Vnpadm(nmd,nmd)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    c1=(reshape((/1,0,0,0,0,0,0,0,0/),shape(c1)))

!vacuum hamiltonian values written directly on the 1/km unit by using PNG(2019values)
    Hvv=reshape((/-0.0879 ,-0.6593 ,-0.7325,-0.6593,-3.4445 ,-0.7325,-3.3150 ,-3.1247,-2.9398/),shape(Hvv)) 


!geometric factor
    Dr=0.5*(((1-sqrt((1-(Eme(3,1)/rtt)**2)))**2)) !Eme(3,1) radius, begining of protoneutronstar
!it is for only after 500 km, equation of motion include only matter and vacuum effect
!matter potential 1/km, bakınız parameter file.


!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (rtt.GT.500) THEN
    do 15 nkk = 1, intme(1,2)

          nkkdz=0+nkk
          u = (intme(1,1)*(nkk+1)*0.5)/intme(1,2)



          dxxfn=(matmul(((Hvv/u)+(yy*0.000271*c1)),psim(nkkdz,:,:) )-matmul( psim(nkkdz,:,:), ((Hvv/u)+(yy*0.000271*c1))) )*(1/ic)
          dxxf(nkkdz,:,:)=dxxfn

!antimatter potential written as negative.
          adxxfn=(matmul(((Hvv/u)-(yy*0.000271*c1)),apsim(nkkdz,:,:) )-matmul( apsim(nkkdz,:,:), ((Hvv/u)-(yy*0.000271*c1)) ))*(1/ic)
          adxxf(nkkdz,:,:)=adxxfn

15        continue
       ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
! integraladm is for adiabatic result. It is used on the odeint part.
          integral = 0.0
          integraladm=0.0
          do 20 nkk=1,intme(1,2)


             nkkdz=0+nkk
             u = (intme(1,1)*(nkk+1)*0.5)/intme(1,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



             nse= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,1)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(1,1));   
             nsm= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,2)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(1,2));
             nst= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,3)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(1,3));
             anse= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,1)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(2,1));
             ansm= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,2)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(2,2));
             anst= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,3)**Eme(3,3)) ))&
                  &* exp(-(Eme(3,3)*u)/Eme(2,3));


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             rhoce=(Eme(3,2)*(1/(Eme(1,1)))*nse)
             rhocm=(Eme(3,2)*(1/(Eme(1,2)))*nsm)
             rhoct=(Eme(3,2)*(1/(Eme(1,3)))*nst)

             arhoce=(Eme(3,2)*(1/(Eme(2,1)))*anse)
             arhocm=(Eme(3,2)*(1/(Eme(2,2)))*ansm)
             arhoct=(Eme(3,2)*(1/(Eme(2,3)))*anst)

             TTN=rhoce+rhocm+rhoct
             ATTN=arhoce+arhocm+arhoct

             radfc=(4.5415d-43*Dr*sqrt(2.0))/(2*pi*(Eme(3,1)**2)) 
             ! it is unitless, self interaction potential unit also 1/km.!don't forget to integral over E      
             !Eme(3,1) radius, begining of protoneutronstar


             valuek(1,1) =( (psim(nkkdz,1,1)*TTN- (apsim(nkkdz,1,1)*ATTN)) )
             valuek(1,2) =( (psim(nkkdz,1,2)*TTN- (apsim(nkkdz,1,2)*ATTN)) )
             valuek(1,3) =( (psim(nkkdz,1,3)*TTN- (apsim(nkkdz,1,3)*ATTN)) )

             valuek(2,1) =( (psim(nkkdz,2,1)*TTN- (apsim(nkkdz,2,1)*ATTN)) )
             valuek(2,2) =( (psim(nkkdz,2,2)*TTN- (apsim(nkkdz,2,2)*ATTN)) )
             valuek(2,3) =( (psim(nkkdz,2,3)*TTN- (apsim(nkkdz,2,3)*ATTN)) )

             valuek(3,1) =( (psim(nkkdz,3,1)*TTN- (apsim(nkkdz,3,1)*ATTN)) )
             valuek(3,2) =( (psim(nkkdz,3,2)*TTN- (apsim(nkkdz,3,2)*ATTN)) )
             valuek(3,3) =( (psim(nkkdz,3,3)*TTN- (apsim(nkkdz,3,3)*ATTN)) )

             if ((nkk.eq.0).or.(nkk.eq.intme(1,2))) then

                integral = integral+valuek
                integraladm = integraladm+valueadm
             else
                integral = integral+(2.0*valuek)
                integraladm = integraladm+(2.0*valueadm)
             end if
20           continue
             hint=intme(1,1)/intme(1,2)

             integral = (hint/2.0)*integral
             integraladm = (hint/2.0)*integraladm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
             Vnp(1,1)=radfc*(integral(1,1))+yy*0.000271
             Vnp(1,2)=radfc*integral(1,2)
             Vnp(1,3)=radfc*integral(1,3)
             Vnp(2,1)=radfc*integral(2,1)
             Vnp(2,2)=radfc*integral(2,2)
             Vnp(2,3)=radfc*integral(2,3)
             Vnp(3,1)=radfc*integral(3,1)
             Vnp(3,2)=radfc*integral(3,2)
             Vnp(3,3)=radfc*integral(3,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1for adiabacity control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             Vnpadm(1,1)=radfc*integraladm(1,1)  !traceless matrix elements (B)
             Vnpadm(2,2)=radfc*integraladm(2,1)  !!traceless matrix elements (B_emuon)
             Vnpadm(3,3)=radfc*integraladm(3,1)

             Vnpadm(1,2)=integraladm(1,1)  !traceless matrix elements (B), derivative
             Vnpadm(1,3)=integraladm(2,1)  !!traceless matrix elements (B_emuon), derivative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             sint=(Vnp)

             do 40 nkk = 1, intme(1,2)

                nkkdz=0+nkk
                u = (intme(1,1)*(nkk+1)*0.5)/intme(1,2)

                dxxfn=(matmul(  ( (Hvv/u)+(sint) ),psim(nkkdz,:,:)  )-matmul(  psim(nkkdz,:,:), ( (Hvv/u) + (sint) ) ))*(1/ic)
                dxxf(nkkdz,:,:)=dxxfn

                adxxfn=(matmul(  ((Hvv/u)-(sint) ),apsim(nkkdz,:,:)  )-matmul(  apsim(nkkdz,:,:), ( (Hvv/u) -(sint) ) ))*(1/ic)

                adxxf(nkkdz,:,:)=adxxfn

40              continue
             END IF
             !write(6,*)t,sint(1,1)
             return
           End Subroutine Derivs


           Subroutine RK4(psiim,apsiim,numb,dxxf,adxxf,t,h,psifm,apsifm,yyh,yyhh,rtt,Tme,Eme,intme,Vnpadm)

             real(kind=kind(8.0D0))  t, h,numb

             real(kind=kind(8.0D0)) th, hh, h6
             double precision :: yyh,yyhh,rtt,rtth
             complex*8,parameter :: ic=( 0.0, 1 )
             integer, parameter :: ned=500, nmd=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
             integer intme(nmd,nmd)

Derivs(yy,rtt,psim,apsim,numb,dxxf,adxxf,Tme,Eme,intme,Vnpadm)

             real(kind=kind(8.0D0)) Tme(nmd,nmd),Eme(nmd,nmd)
             complex(kind=kind(8.0D0)) dxxf(ned,nmd,nmd),dxxfm(ned,nmd,nmd),dxfx(ned,nmd,nmd)
             complex(kind=kind(8.0D0)) adxxf(ned,nmd,nmd),adxxfm(ned,nmd,nmd),adxfx(ned,nmd,nmd)
             complex(kind=kind(8.0D0))psim(ned,nmd,nmd) ,psiim(ned,nmd,nmd), psifm(ned,nmd,nmd)   
             complex(kind=kind(8.0D0))apsim(ned,nmd,nmd) ,apsiim(ned,nmd,nmd), apsifm(ned,nmd,nmd),Vnpadm(nmd,nmd)   
!Given values for the variables psiim(ned,nmd,nmd) and their derivatives dxxf(1:n) known at rtt , use
!the fourth-order Runge-Kutta method to advance the solution over an interval h and return
!the incremented variables as psifm(ned,nmd,nmd) , which need not be a distinct array from psiim . The
!user supplies the subroutine derivs, which returns derivatives dxxf at rtt .

             hh = 0.5d0 * h
             h6 =h/6.d0
             th = t + hh

             rtth=rtt+hh

             psim = psiim+ hh*dxxf
             apsim = apsiim+ hh*adxxf

             call Derivs(yyhh,rtth,psim,apsim,numb,dxfx,adxfx,Tme,Eme,intme,Vnpadm)
             psim= psiim + hh * dxfx
             apsim = apsiim + hh * adxfx
             call Derivs(yyhh,rtth,psim,apsim,numb,dxxfm,adxxfm,Tme,Eme,intme,Vnpadm)
             psim = psiim+ h*dxxfm
             apsim = apsiim+ h*adxxfm
             dxxfm= dxfx + dxxfm
             adxxfm= adxfx + adxxfm
             call Derivs(yyh,rtt+h,psim,apsim,numb,dxfx,adxfx,Tme,Eme,intme,Vnpadm)



             psifm =( psiim+h6*(dxxf+dxfx+2.d0*dxxfm))
             apsifm=( apsiim+h6*(adxxf+adxfx+2.d0*adxxfm))
           End Subroutine RK4


           Subroutine rkqc(psiim,apsiim,numb,yy,yyh,yyhh,rtt,dxxf,adxxf, n, t, htry, eps,psiscal,apsiscal, hdid, hnext,&
                &Tme,Eme,intme,Vnpadm)
!          Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracyand  adjust  stepsize.  Input  are  the  dependent
!          variable  vectory(1:n)and  its  derivativedydx(1:n)at the starting value of the independent variablex.  Also input are the stepsizeto  be 
!         attemptedhtry,  the  required  accuracyeps,  and  the  vectoryscal(1:n)againstwhich the error is scaled.  On output,yandxare replaced by their new
! values,hdidis thestepsize that was actually accomplished, andhnextis the estimated next stepsize.derivsis the user-supplied subroutine that computes the
! right-hand side derivative

             real(kind=kind(8.0D0))  t, htry, eps, hdid, hnext,numb

             real(kind=kind(8.0D0)) tsav,hh,h,temp,errmax

             real(kind=kind(8.0D0)) pgron, pshrnk,fcor,un,safety,errcon,tiny   

             integer, parameter :: ned=500, nmd=3
             integer intme(nmd,nmd),n


             real(kind=kind(8.0D0))Tme(nmd,nmd),Eme(nmd,nmd)
             complex(kind=kind(8.0D0))dxxf(ned,nmd,nmd) ,dxsav(ned,nmd,nmd),&
                  & adxxf(ned,nmd,nmd) ,adxsav(ned,nmd,nmd)

             complex(kind=kind(8.0D0))psiim(ned,nmd,nmd),apsiim(ned,nmd,nmd)
             complex(kind=kind(8.0D0))psiscal(ned,nmd,nmd),psisav(ned,nmd,nmd),&
                  &apsisav(ned,nmd,nmd),apsiscal(ned,nmd,nmd),Vnpadm(nmd,nmd)
             complex(kind=kind(8.0D0))psitemp(ned,nmd,nmd),apsitemp(ned,nmd,nmd)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
             integer :: i,nk,nkk,nkkd,nin,end_val
             double precision :: yy,rtt,yyh,yyhh

             pgron=-0.20d0
             pshrnk=-0.25d0
             fcor=0.06666666d0   !1/15
             un = 1.d0;
             safety=0.9d0
             errcon=6.D-4
             tiny= 1.D-20   

             tsav= t          !Save begin time



             psisav = psiim
             apsisav = apsiim

             dxsav= dxxf
             adxsav= adxxf




             h = htry         !define increment for a try value
1            hh = 0.5d0*h     !take 2 half time steps





             call rk4(psisav,apsisav,numb,dxsav,adxsav,tsav,hh,psitemp,apsitemp,yyh,yyhh,rtt,Tme,Eme,intme,Vnpadm)

             t = tsav + hh

             call Derivs(yy,rtt,psitemp,apsitemp,numb,dxxf,adxxf,Tme,Eme,intme,Vnpadm)
             call rk4(psitemp,apsitemp,numb,dxxf,adxxf,t,hh,psiim,apsiim,yyh,yyhh,rtt,Tme,Eme,intme,Vnpadm)
             t = tsav + h

             !         write(6,103) t, ytemp2(1),ytemp2(2) 
             if (t == tsav) then
                print *,' Pause in RKQC subroutine'
                print *,' Increment too small of independant variable'
                !           pause ' Press any key to continue...'
             end if

             call rk4(psisav,apsisav,numb,dxsav,adxsav,tsav,h,psitemp,apsitemp,yyh,yyhh,rtt,Tme,Eme,intme,Vnpadm)
             !        write(6,103) tsav, ytemp2(1),ytemp2(2)
             !          print *,tsav, ytemp2(1),ytemp2(2)
             errmax = 0.d0    !Evaluate error
             temp = 0.d0

             do i = 1,n
                do nkk = 1, intme(1,2) 

                   nkkd=0+nkk

                   psitemp( nkkd,i,1) = psiim( nkkd,i,1) - psitemp( nkkd,i,1)     !psitemp = estimated error

                   IF (abs(real(psiscal(nkkd,i,1)))>tiny)  temp = abs(psitemp( nkkd,i,1)/psiscal(nkkd,i,1)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do

             end do

             do i = 1,n
                do nkk = 1, intme(1,2) 

                   nkkd=0+nkk
                   apsitemp( nkkd,i,1) = apsiim( nkkd,i,1) - apsitemp( nkkd,i,1)     !psitemp = estimated error
                   IF (abs(real(apsiscal(nkkd,i,1)))>tiny)  temp = abs(apsitemp( nkkd,i,1)/apsiscal(nkkd,i,1)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do
             end do

             do i = 1,n
                do nkk = 1, intme(1,2) 

                   nkkd=0+nkk
                   psitemp( nkkd,i,2) = psiim( nkkd,i,2) - psitemp( nkkd,i,2)    !psitemp = estimated error

                   IF (abs(real(psiscal(nkkd,i,2)))>tiny)  temp = abs(psitemp( nkkd,i,2)/psiscal(nkkd,i,2)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do
             end do

             do i = 1,n
                do nkk = 1, intme(1,2) 

                   nkkd=0+nkk
                   apsitemp( nkkd,i,2) = apsiim( nkkd,i,2) - apsitemp( nkkd,i,2)    !psitemp = estimated error 
                   IF (abs(real(apsiscal(nkkd,i,2)))>tiny)  temp = abs(apsitemp( nkkd,i,2)/apsiscal(nkkd,i,2)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do
             end do

             do i = 1,n
                do nkk = 1, intme(1,2) 

                   nkkd=0+nkk
                   psitemp( nkkd,i,3) = psiim( nkkd,i,3) - psitemp( nkkd,i,3)    !psitemp = estimated error

                   IF (abs(real(psiscal(nkkd,i,3)))>tiny)  temp = abs(psitemp( nkkd,i,3)/psiscal(nkkd,i,3)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do
             end do

             do nkk = 1, intme(1,2)

                nkkd=0+nkk
                do i = 1,n
                   apsitemp( nkkd,i,3) = apsiim( nkkd,i,3) - apsitemp( nkkd,i,3)    !psitemp = estimated error

                   IF (abs(real(apsiscal(nkkd,i,3)))>tiny)  temp = abs(apsitemp( nkkd,i,3)/apsiscal(nkkd,i,3)) !there is a real
                   IF (errmax < temp)  errmax = temp
                end do
             end do

             errmax = errmax/eps             !real error / requirement
             if (errmax > un) then           !Error too big, reduce h
                h = safety*h*exp(pshrnk*dlog(errmax))
                GOTO 1                       !start again
             else                            !the step has been a success
                hdid = h                     !Calculate next time step
                IF (errmax > errcon) THEN
                   hnext = safety*h*dexp(pgron*dlog(errmax))
                ELSE
                   hnext = 4.d0*h
                end if
             end if

             psiim = psiim + psitemp*fcor
             apsiim = apsiim + apsitemp*fcor




             return
             !103      format(3(1pe12.3))
           End Subroutine rkqc

           Subroutine odeint(psistart,apsistart, nvar,numb,xxx,yyy,nmax, t1, t2, eps, h1, hmin, nok, nbad,&
                &Tme,Eme,intme,Vnpadm)
!Runge-Kutta driver with adaptive stepsize control.  Integrate the starting valuesystart(1:nvar)fromx1tox2with accuracyeps, storing intermediate results in !the common block/path/.h1should  be set as a guessed  first stepsize,hminas the minimum allowed stepsize (canbe  zero).  On  outputnokandnbadare  the  number  of  good  and  bad  (but  retried and
!ixed) steps taken, andystartis replaced by values at the end of the integration interval.derivsis the user-supplied subroutine for calculating the right-hand side derivative, while rkqc is the name of the stepper routine to be used./path/contains its own informationabout  how often an intermediate  value is to be stored.

             integer, parameter :: ned=500, nmd=3
             integer nstp,maxstp
             real(kind=kind(8.0D0)) t1, t2, eps, h1, hmin,numb,valuwrit,nok,nbad
             integer intme(nmd,nmd)

             real*8 two,zero,tiny,tsav,t,hnext,hdid,h,dtsave,R,F2,g,u,r_cm
             !complex*8 yscal(1:nvar),y(1:nvar),dydt(50,50),ystart(*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111


             complex(kind=kind(8.0D0)) dxxf(ned,nmd,nmd),adxxf(ned,nmd,nmd),UM(nmd,nmd),MC(nmd,nmd),AMC(nmd,nmd)
             complex(kind=kind(8.0D0)) psiim(ned,nmd,nmd),apsiim(ned,nmd,nmd),masspsiim (ned,nmd,nmd),&
                  &amasspsiim (ned,nmd,nmd),Vnpadm(nmd,nmd)               
             complex(kind=kind(8.0D0))psiscal(ned,nmd,nmd),psistart(ned,nmd,nmd),newmatr(nmd,nmd)  
             complex(kind=kind(8.0D0))apsiscal(ned,nmd,nmd),apsistart(ned,nmd,nmd) ,anewmatr(nmd,nmd)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
             integer :: k, nmax,jlo,imax,nvar,nkkdz,nkkz,nin,end_val
             double precision :: yy,yyh,yyhh,dyy,rtt,dyyh,dyyhh
             !double precision, dimension(5,200) ::  datei
             double precision, dimension(10000000) :: xxx,yyy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
             real*8 nse,nsm,nst,anse,ansm,anst, Te,Tm,Tt,aTe,aTm,aTt
             real*8 aEe,aEm,aEt,Ee,Em,Et,L_g,TTN,ATTN
             real*8 rhoce,rhocm,rhoct,arhoce,arhocm,arhoct
             real(kind=kind(8.0D0)) Tme(nmd,nmd),Eme(nmd,nmd)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

             real, parameter :: pi=3.141592654
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
             real*8 omega,delta_m_sqr31,theta13,omegcos,omegsinsqr,omegasin
             real*8 eigenval,theta,tht,mass1,mass2,sintheta,n_edde,ratioelcden
             real*8 omegam,omegcosm,omegsinsqrm,omegasinm,n_ed,n_eres

             real*8 m1,m2,mattm1,mattm2,itm,stm,potderiv,adiabacity,radfch,Drh
             real*8 selfederiv,selfemuderiv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             maxstp= 1000000000
             two   =    2
             zero  =    0
             tiny  = 1.d-20
             dtsave = 0.d0       !not used here

             t= t1
             IF (t2 > t1) THEN 
                h =   abs(h1)
             ELSE 
                h = - abs(h1)
             END IF
             nok = 0
             nbad = 0

             psiim=psistart
             apsiim=apsistart

             tsav = t - dtsave*two

             do nstp = 1, maxstp
                r_cm=t
                call hunt(xxx,nmax,r_cm,jlo)



                k=min(max(jlo-1/2,1),nmax -1) !for 2-point interpolation

                call polint(xxx(k),yyy(k),2,r_cm,r_cm+h,r_cm+0.5*h,yy,yyh,yyhh,dyy,dyyh,dyyhh)



                rtt=r_cm


                !print*,rtt,h




                call Derivs(yy,rtt,psiim,apsiim,numb,dxxf,adxxf,Tme,Eme,intme,Vnpadm)
                !print*,r_cm,yy

                psiscal = abs(psiim)+abs(dxxf*h)
                apsiscal = abs(apsiim)+abs(adxxf*h)


                valuwrit=t-numb

                IF (valuwrit .GT.200) THEN

                   numb=rtt

                   do nkkz = 1, intme(1,2)

                      nkkdz=0+nkkz

                      u = (intme(1,1)*(nkkz+1)*0.5)/intme(1,2)

                      nse= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,1)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(1,1));   !!!!loook at!!!!!!!
                      nsm= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,2)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(1,2));
                      nst= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(1,3)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(1,3));
                      anse= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,1)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(2,1));
                      ansm= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,2)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(2,2));
                      anst= (((Eme(3,3)**Eme(3,3))/6)*((u**(Eme(3,3)-1))/(Eme(2,3)**Eme(3,3)) ))&
                           &* exp(-(Eme(3,3)*u)/Eme(2,3));

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      rhoce=(Eme(3,2)*(1/(Eme(1,1)))*nse)
                      rhocm=(Eme(3,2)*(1/(Eme(1,2)))*nsm)
                      rhoct=(Eme(3,2)*(1/(Eme(1,3)))*nst)


                      arhoce=(Eme(3,2)*(1/(Eme(2,1)))*anse)
                      arhocm=(Eme(3,2)*(1/(Eme(2,2)))*ansm)
                      arhoct=(Eme(3,2)*(1/(Eme(2,3)))*anst)

                      TTN=rhoce+rhocm+rhoct

                      ATTN=arhoce+arhocm+arhoct


                      !Drh=0.5*(((1-sqrt((1-(18/(rtt+h))**2)))**2))
                      !radfch=(4.5415d-43*Drh*sqrt(2.0))/(2*pi*(Eme(3,1)**2))    !Eme(3,1) radius, begining of protoneutronstar

                      !itm=(delta_m_sqr31*cos(2*theta13)+((2.0*u)*( yy*0.000271)+(real(Vnpadm(1,1)))))  !potential negative

                      !theta=(acos((delta_m_sqr31*cos(0.2940)+(2*u*(yy*0.000271)))/eigenval)*0.5) !potential negative

                      !stm=(  delta_m_sqr31*sin(2*theta13) + (4*u*real(Vnpadm(2,2)))     )
                      !eigenval=(sqrt( (itm)**2+(stm)**2) )
                      !theta=(acos((delta_m_sqr31*cos(2*theta13)+(2*u*(yy*0.000271-2*real(Vnpadm(1,1)))))/eigenval)*0.5) 

                      !mass1=real(((cos(theta)**2)*apsiim(nkkdz,1,1))+( (sin(theta)**2)*2*apsiim(nkkdz,2,2))&
                      !&-(apsiim(nkkdz,1,3))*cos(theta)*sin(theta)-(apsiim(nkkdz,3,1)*sin(theta)*cos(theta)))

                      !mass2=real(((sin(theta)**2)*apsiim(nkkdz,1,1))+( (cos(theta)**2)*2*apsiim(nkkdz,2,2))&
                      !&+(apsiim(nkkdz,1,3))*cos(theta)*sin(theta)+(apsiim(nkkdz,3,1)*sin(theta)*cos(theta)))

                      !potderiv= abs(((yyh*0.000271)-(yy*0.000271))/h)

                      !selfederiv=abs ((real(Vnpadm(1,2))*radfch-real(Vnpadm(1,1)))/h)    !dioganal element
                      !selfemuderiv=abs (real((Vnpadm(1,3))*radfch-real(Vnpadm(2,2)))/h)  !nondiagonal element

                      !adiabacity= ( (eigenval)**2)/( 4*(u**2)*(  sin(2*theta)*( potderiv+(2*selfederiv) )&
                      !&+ ( 2*cos(2*theta)*(2*selfemuderiv) )))

                      !n_ed= (-yy*(10d+15)*0.5)/(1.67*(10d-24))


                      !n_eres=((delta_m_sqr31)*cos(2*theta13))/(2*sqrt(2.0)*(4.5415d-43)*u)

                      !m1=(0.15*1e-15)*(1/(1.97*1e-16))  !square of m1

                      !m2= (2.50*1e-15)*(1/(1.97*1e-16))    !square of m2

                      !mattm1=(((0.5)*( (m1)+(m2)+((-yy*0.000271)*u*2)+eigenval))*(1.97*(1d-4)))
                      !mattm2=(((0.5)*((m1)+(m2)+((-yy*0.000271)*u*2)-eigenval))*(1.97*(1d-4)))



                      WRITE (UNIT=14, FMT=*)u,t,real(psiim(nkkdz,1,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,1,2))*(TTN*Tme(3,1)),&
                           &real(psiim(nkkdz,1,3))*(TTN*Tme(3,1)),real(psiim(nkkdz,2,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,2,2))*(TTN*Tme(3,1)),&
                           &real(psiim(nkkdz,2,3))*(TTN*Tme(3,1)),real(psiim(nkkdz,3,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,3,2))*(TTN*Tme(3,1)),&
                           &real(psiim(nkkdz,3,3))*(TTN*Tme(3,1))

                      WRITE (UNIT=15, FMT=*)u,t,real(apsiim(nkkdz,1,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,1,2))*(ATTN*Tme(3,1)),&
                           &real(apsiim(nkkdz,1,3))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,2,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,2,2))*(ATTN*Tme(3,1)),&
                           &real(apsiim(nkkdz,2,3))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,3,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,3,2))*(ATTN*Tme(3,1)),&
                           &real(apsiim(nkkdz,3,3))*(ATTN*Tme(3,1))

                      WRITE (UNIT=16, FMT=*)u,t,real(psiim(nkkdz,1,1)),real(psiim(nkkdz,1,2)),&
                           &real(psiim(nkkdz,1,3)),real(psiim(nkkdz,2,1)),real(psiim(nkkdz,2,2)),&
                           &real(psiim(nkkdz,2,3)),real(psiim(nkkdz,3,1)),real(psiim(nkkdz,3,2)),&
                           &real(psiim(nkkdz,3,3))

                      WRITE (UNIT=17, FMT=*)u,t,real(apsiim(nkkdz,1,1)),real(apsiim(nkkdz,1,2)),&
                           &real(apsiim(nkkdz,1,3)),real(apsiim(nkkdz,2,1)),real(apsiim(nkkdz,2,2)),&
                           &real(apsiim(nkkdz,2,3)),real(apsiim(nkkdz,3,1)),real(apsiim(nkkdz,3,2)),&
                           &real(apsiim(nkkdz,3,3))


                      !WRITE (UNIT=18, FMT=*)u,t,eigenval,theta,sin(2.0*theta),&
                      !&mass1,mass2,sin(2*theta),real(psiim(nkkdz,1,1)),real(psiim(nkkdz,2,2)),&
                      !&real(apsiim(nkkdz,1,1)),real(apsiim(nkkdz,2,2)),mattm1,mattm2,(n_ed/n_eres),&
                      !&adiabacity,real(Vnpadm(1,1)),real(Vnpadm(2,2))
                      !,yy,h,eigenval,theta,u,potderiv
                      !((eigenval)**2)/(4*(u**2)*sin(2*theta)*potderiv)

                      !WRITE (UNIT=19, FMT=*)r_cm,yy,r_cm+h,yyh,r_cm+0.5*h,yyhh

                   end do

                endif


                IF (((t+h-t2)*(t+h-t1)) > zero)   h = t2 - t
                call rkqc(psiim,apsiim,numb,yy,yyh,yyhh,rtt,dxxf,adxxf,nvar,t,h,eps,psiscal,apsiscal, hdid,hnext,&
                     &Tme,Eme,intme,Vnpadm)

                IF (hdid == h ) THEN 
                   nok = nok + 1
                ELSE 
                   nbad = nbad + 1
                END IF
                IF (((t-t2)*(t2-t1)) >= zero) THEN

                   psistart=psiim
                   apsistart=apsiim


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   GOTO 99  !it is over
                END IF
                IF (abs(hnext) < hmin) THEN
                   nok = -1   !error flag
                   !           pause ' Time step too small!'
                   GOTO 99
                END IF
                h = hnext
             end do

             !	 pause ' Pause in subroutine ODEINT - too many time steps!'
99           h1 = h

             WRITE (UNIT=18, FMT=*)u,t,real(psiim(nkkdz,1,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,1,2))*(TTN*Tme(3,1)),&
                  &real(psiim(nkkdz,1,3))*(TTN*Tme(3,1)),real(psiim(nkkdz,2,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,2,2))*(TTN*Tme(3,1)),&
                  &real(psiim(nkkdz,2,3))*(TTN*Tme(3,1)),real(psiim(nkkdz,3,1))*(TTN*Tme(3,1)),real(psiim(nkkdz,3,2))*(TTN*Tme(3,1)),&
                  &real(psiim(nkkdz,3,3))*(TTN*Tme(3,1))

             WRITE (UNIT=19, FMT=*)u,t,real(apsiim(nkkdz,1,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,1,2))*(ATTN*Tme(3,1)),&
                  &real(apsiim(nkkdz,1,3))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,2,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,2,2))*(ATTN*Tme(3,1)),&
                  &real(apsiim(nkkdz,2,3))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,3,1))*(ATTN*Tme(3,1)),real(apsiim(nkkdz,3,2))*(ATTN*Tme(3,1)),&
                  &real(apsiim(nkkdz,3,3))*(ATTN*Tme(3,1))



             !return
           End Subroutine odeint

           subroutine hunt(xx,n,x,jlo)
             !(xxx,nmax,r_cm(ij),r_cm(ij)+h,jlo)
             !----------linkages.
             !     called by - [subroutine] start, trajectory, rate_reaclib, rate5, rate6
             !     calls     - none

             !     This subroutine performs hunting of an appropriate data value.
             !     Based on Numerical Recipes, Cambridge University Press; 2nd edition (1992).

             implicit none
             integer jlo,n
             double precision x,xx(n)
             integer inc,jhi,jm
             logical ascnd

             ascnd=xx(n).ge.xx(1)
             if(jlo.le.0.or.jlo.gt.n)then
                jlo=0
                jhi=n+1
                go to 3
             end if

             inc=1
             if(x.ge.xx(jlo).eqv.ascnd)then
1               jhi=jlo+inc
                if(jhi.gt.n)then
                   jhi=n+1
                else if(x.ge.xx(jhi).eqv.ascnd)then
                   jlo=jhi
                   inc=inc+inc
                   go to 1
                end if
             else
                jhi=jlo
2               jlo=jhi-inc
                if(jlo.lt.1)then
                   jlo=0
                else if(x.lt.xx(jlo).eqv.ascnd)then
                   jhi=jlo
                   inc=inc+inc
                   go to 2
                end if
             end if

3            if(jhi-jlo.eq.1)then
                if(x.eq.xx(n)) jlo=n-1
                if(x.eq.xx(1)) jlo=1
                return
             end if

             jm=(jhi+jlo)/2
             if(x.ge.xx(jm).eqv.ascnd)then
                jlo=jm
             else
                jhi=jm
             end if
             go to 3
             !return      
           end subroutine hunt

           !call polint(xxx(k),xxx(k)+0.5d0*h,yyy(k),2,r_cm(ij),yy,yyh,dyy,dyyh)
           !===========================================================================
           subroutine polint(xa,ya,n,x,xh,xhh,y,yh,yhh,dy,dyh,dyhh)
             !              (xxx(k),yyy(k),2,r_cm(ij),r_cm(ij)+h,r_cm(ij)+0.5*h,yy,yyh,yyhh,dyy,dyyh,dyyhh)
             !----------linkages.
             !     called by - [subroutine] start, rate_reaclib, trajectory, rate5, rate6
             !     calls     - none
             !     This subroutine performs polynomial interpolation and extrapolation in normal scale space.
             !     Based on Numerical Recipes, Cambridge University Press; 2nd edition (1992).

             implicit none
             integer n,nmax
             double precision dy,x,y,xa(n),ya(n),xh,yh,dyh,xhh,yhh,dyhh
             parameter (nmax=10)
             integer i,m,ns
             double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
             double precision denh,difh,difth,hoh,hph,wh,ch(nmax),dh(nmax)
             double precision denhh,difhh,difthh,hohh,hphh,whh,chh(nmax),dhh(nmax)
             !-----interpolation or extrapolation
             ns=1

             dif=abs(x-xa(1))
             difh=abs(xh-xa(1))
             difhh=abs(xhh-xa(1))
             do i=1,n

                dift=abs(x-xa(i))
                difth=abs(xh-xa(i))
                difthh=abs(xhh-xa(i))

                if(dift.lt.dif)then
                   ns=i

                   dif=dift
                   difh=difth
                   difhh=difthh

                end if
                c(i)=ya(i)
                d(i)=ya(i)

                ch(i)=ya(i)
                dh(i)=ya(i)

                chh(i)=ya(i)
                dhh(i)=ya(i)
             end do

             y=ya(ns)

             yh=ya(ns)
             yhh=ya(ns)
             ns=ns-1
             do m=1,n-1
                do i=1,n-m

                   ho=xa(i)-x
                   hp=xa(i+m)-x

                   hoh=xa(i)-xh
                   hph=xa(i+m)-xh

                   hohh=xa(i)-xhh
                   hphh=xa(i+m)-xhh


                   w=c(i+1)-d(i)

                   wh=ch(i+1)-dh(i)
                   whh=chh(i+1)-dhh(i)

                   den=ho-hp
                   denh=hoh-hph
                   denhh=hohh-hphh

                   if(den.eq.0.)then
                      write(*,*) i,m,xa(i),xa(i+m),x
                      write(*,*) 'pause:  failure in polint'
                      read(*,*)
                   end if

                   den=w/den
                   denh=wh/denh
                   denhh=whh/denhh


                   d(i)=hp*den
                   c(i)=ho*den

                   dh(i)=hph*denh
                   ch(i)=hoh*denh

                   dhh(i)=hphh*denhh
                   chh(i)=hohh*denhh

                end do

                if(2*ns.lt.n-m)then
                   dy=c(ns+1)

                   dyh=ch(ns+1)
                   dyhh=chh(ns+1)

                else
                   dy=d(ns)
                   dyh=dh(ns)
                   dyhh=dhh(ns)

                   ns=ns-1
                end if
                y=y+dy
                yh=yh+dyh
                yhh=yhh+dyhh

             end do

             !  return
           end subroutine polint






