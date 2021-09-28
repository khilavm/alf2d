! Gridsetting subroutine. 2-D Yee grid to solve Maxwell's equations.

! collfreq.txt has the ion-neutral collision frequencies, plasmafreq1 is for the plasma frequencies of electrons and ions, gyrofreq has the gyrofrequencies, wavespeed has the final speed of the travelling wave and pedersen,parallel_conds has the conductivities to be used in solving the wave equations.

! Speed of light is in km/s.

! The data for ions is in terms of percentages, so the electron density has been used to convert the percentages to densities. This is easy to do as all ions are singly ionised.

! Wave speed and Pedersen conductivity calculated from formula in RL '99, parallel conductivity from that given in his space physics notes.


program grid
implicit none

integer :: i,k,j,m
real(8),allocatable,dimension(:) :: height,collfreq,etemp,edens,odens,hdens,hedens,o2dens,nodens,ndens
real(8),allocatable,dimension(:) :: nuo,nun2,nuo2,nuhe,nuar,nuh,nun,b,ne,beta,pedrev,rveden,abc,bbb
real(8),allocatable,dimension(:) :: pfe,pfo,pfh,pfhe,pfo2,pfno,pfn,gye,gyo,gyh,gyhe,gyo2,gyno,gyn,alpha,speed
real(8),allocatable,dimension(:) :: cfoo,cfoh,cfohe,cfoo2,cfono,cfon,cfo2o,cfo2h,cfo2he,cfo2o2,cfo2no,cfo2n,cfei
real(8),allocatable,dimension(:) :: cfn2o,cfn2h,cfn2he,cfn2o2,cfn2no,cfn2n,sigped,sigpar,rvcfei,epspar
real(8),allocatable,dimension(:) :: cfheo,cfheh,cfhehe,cfheo2,cfheno,cfhen,cfaro,cfarh,cfarhe,cfaro2,cfarno,cfarn
real(8),allocatable,dimension(:) :: cfho,cfhh,cfhhe,cfho2,cfhno,cfhn,cfno,cfnh,cfnhe,cfno2,cfnno,cfnn
real(8),parameter :: echarge=1.6E-19,eps=8.85E-12,mprot=1.67E-27,melec=9.1E-31,c=3.0E5,j0=1.0E-6
real(8),parameter :: const=(echarge**2)/mprot/eps,const1=(echarge**2)*(1.0E-18)/(mprot)**2
real(8),parameter :: const2=(echarge**2)*(1.0E-18)/(melec)**2,const3=2.6E-9,c4=(echarge**2)/mprot,c5=(echarge**2)/melec

real(8),dimension(950) :: edens1,odens1,hdens1,ndens1,hedens1,o2dens1,nodens1
real(8),dimension(950) :: z,b1,cfoo1,cfoh1,cfohe1,cfoo21,cfono1,cfon1,cfo2o1,cfo2h1,cfo2he1,cfo2o21,cfo2no1,cfo2n1
real(8),dimension(950) :: cfn2o1,cfn2h1,cfn2he1,cfn2o21,cfn2no1,cfn2n1,cfei1,sigped1
real(8),dimension(950) :: cfheo1,cfheh1,cfhehe1,cfheo21,cfheno1,cfhen1,cfaro1,cfarh1,cfarhe1,cfaro21,cfarno1,cfarn1
real(8),dimension(950) :: cfho1,cfhh1,cfhhe1,cfho21,cfhno1,cfhn1,cfno1,cfnh1,cfnhe1,cfno21,cfnno1,cfnn1
real(8),dimension(996) :: b2,hx,hy,edens2,ne1,odens2,hdens2,hedens2,o2dens2,nodens2,ndens2,dentot1
real(8),dimension(996) :: cfoo_2,cfoh_2,cfohe_2,cfoo2_2,cfono_2,cfon_2
real(8),dimension(996) :: cfo2o_2,cfo2h_2,cfo2he_2,cfo2o2_2,cfo2no_2,cfo2n_2
real(8),dimension(996) :: cfn2o_2,cfn2h_2,cfn2he_2,cfn2o2_2,cfn2no_2,cfn2n_2
real(8),dimension(996) :: cfheo_2,cfheh_2,cfhehe_2,cfheo2_2,cfheno_2,cfhen_2
real(8),dimension(996) :: cfaro_2,cfarh_2,cfarhe_2,cfaro2_2,cfarno_2,cfarn_2
real(8),dimension(996) :: cfho_2,cfhh_2,cfhhe_2,cfho2_2,cfhno_2,cfhn_2
real(8),dimension(996) :: cfno_2,cfnh_2,cfnhe_2,cfno2_2,cfnno_2,cfnn_2
real(8),dimension(996) :: pfe1,pfo1,pfh1,pfhe1,pfo21,pfno1,pfn1,gye1,gyo1,gyh1,gyhe1,gyo21,gyno1,gyn1
real(8),dimension(996) :: alpha1,speed1,va1,sigped2,cfei_2,sigpar1,epspar1,beta1,pedrev1,rveden1,rvcfei1,rva1
real(8),dimension(996) :: aaa1,bb1,bbb1,abc1,cc1,bbb2,cour1,rvepspar1,sigped_1
real(8),dimension(12,996) :: Jz2,Ex2,By2,Ez2

real(8),allocatable,dimension(:) :: cour,aaa,bb,cc,ex1,x,dentot,va,rva,jdrive,by1
real(8),allocatable,dimension(:,:) :: Jz,Ex,By,Ez
real(8),parameter :: dt=1.0E-3,dx=5,dz=20,mu=1.26E-6,tramp=200,ex0=3.0E-3,jz0=1.0E-6,re=6400
real(8) :: t,aa,ab,hintpedz,hintpeds,dump1,dump2,dump3,dump4
integer(8) :: n

allocate(height(46),etemp(46),edens(46),collfreq(46),odens(46),hdens(46),hedens(46),o2dens(46),nodens(46),ndens(46))
allocate(nuo(46),nun2(46),nuo2(46),nuhe(46),nuar(46),nuh(46),nun(46),b(46),beta(46),rvcfei(46))
allocate(cfn2o(46),cfn2h(46),cfn2he(46),cfn2o2(46),cfn2no(46),cfn2n(46),cfoo(46),cfoh(46),cfohe(46),cfoo2(46))
allocate(cfono(46),cfon(46),cfo2o(46),cfo2h(46),cfo2he(46),cfo2o2(46),cfo2no(46),cfo2n(46),cfei(46),cfnn(46))
allocate(cfheo(46),cfheh(46),cfhehe(46),cfheo2(46),cfheno(46),cfhen(46),cfaro(46),cfarh(46),cfarhe(46),cfaro2(46),cfarno(46))
allocate(cfarn(46),cfho(46),cfhh(46),cfhhe(46),cfho2(46),cfhno(46),cfhn(46),cfno(46),cfnh(46),cfnhe(46),cfno2(46),cfnno(46))
allocate(pfe(46),pfo(46),pfh(46),pfhe(46),pfo2(46),pfno(46),pfn(46),sigped(46),sigpar(46),pedrev(46),rveden(46))
allocate(gye(46),gyo(46),gyh(46),gyhe(46),gyo2(46),gyno(46),gyn(46),alpha(46),speed(46),ne(46))
allocate(dentot(46),va(46),rva(46))
allocate(cour(46),aaa(46),bb(46),cc(46),ex1(12),x(12))
allocate(epspar(46),abc(46),bbb(46),jdrive(12),by1(12))
allocate(By(12,996),Ex(12,996),Jz(12,996),Ez(12,996))

open (10,file="Etemp_1.txt")
open (11,file="Edens,O,H,He,O2,NO,N_1.txt")
open (12,file='O,N2,O2,He,Ar,H,N_1.txt')
open (13,file='B field_1.txt')
!open (1,file="collfreq_2.txt",status='new')
!open (2,file="plasmafreq_2.txt",status='new')
!open (3,file='gyrofreq_2.txt',status='new')
!open (4,file='wavespeed_2.txt',status='new')
!open (5,file='pedersen,parallel_conds_2.txt',status='new')
!open (6,file='sigpeds.dat',status='new')
open (90,file='fields2_100.csv',status='new')
open (91,file='fields2_200.csv',status='new')
open (92,file='fields2_300.csv',status='new')
open (93,file='fields2_400.csv',status='new')
open (94,file='fields2_500.csv',status='new')
open (95,file='fields2_600.csv',status='new')
open (96,file='fields2_700.csv',status='new')
open (97,file='fields2_800.csv',status='new')
open (98,file='fields2_900.csv',status='new')
open (99,file='fields2_980.csv',status='new')
open (80,file='fields5_100.csv',status='new')
open (81,file='fields5_200.csv',status='new')
open (82,file='fields5_300.csv',status='new')
open (83,file='fields5_400.csv',status='new')
open (84,file='fields5_500.csv',status='new')
open (85,file='fields5_600.csv',status='new')
open (86,file='fields5_700.csv',status='new')
open (87,file='fields5_800.csv',status='new')
open (88,file='fields5_900.csv',status='new')
open (89,file='fields5_980.csv',status='new')
open (70,file='fields8_100.csv',status='new')
open (71,file='fields8_200.csv',status='new')
open (72,file='fields8_300.csv',status='new')
open (73,file='fields8_400.csv',status='new')
open (74,file='fields8_500.csv',status='new')
open (75,file='fields8_600.csv',status='new')
open (76,file='fields8_700.csv',status='new')
open (77,file='fields8_800.csv',status='new')
open (78,file='fields8_900.csv',status='new')
open (79,file='fields8_980.csv',status='new')
open (60,file='fields11_100.csv',status='new')
open (61,file='fields11_200.csv',status='new')
open (62,file='fields11_300.csv',status='new')
open (63,file='fields11_400.csv',status='new')
open (64,file='fields11_500.csv',status='new')
open (65,file='fields11_600.csv',status='new')
open (66,file='fields11_700.csv',status='new')
open (67,file='fields11_800.csv',status='new')
open (68,file='fields11_900.csv',status='new')
open (69,file='fields11_980.csv',status='new')

do i=1,46
	read (10,fmt='(2g11.3)') height(i),etemp(i)
	read (12,fmt='(7g12.3)') nuo(i),nun2(i),nuo2(i),nuhe(i),nuar(i),nuh(i),nun(i) ! neutral densities (cm^-3)
	read (11,fmt='(7g13.5)') edens(i),odens(i),hdens(i),hedens(i),o2dens(i),nodens(i),ndens(i) ! ion densities as percentage of eletron density. electron density is in m^-3.
	read (13,fmt='(1F7.1)') b(i)
	ne(i)=edens(i)*(1.0E-6) ! converting to cm^-3.
	dentot(i)=(melec*edens(i))+(mprot*16*odens(i)*edens(i)/100)+(mprot*hdens(i)*edens(i)/100)+ &
			(mprot*4*hedens(i)*edens(i)/100)+(mprot*32*o2dens(i)*edens(i)/100)+ &
			(mprot*30*nodens(i)*edens(i)/100)+(mprot*14*ndens(i)*edens(i)/100) ! in SI units
	!! Calculation of ion-neutral collision frequencies^2 follows. Formula requires densities in cm^-3.
	cfoo(i)=(const3*(nuo(i)+(ne(i)*odens(i)/100))/4)**2
	cfoh(i)=(const3*(nuo(i)+(ne(i)*hdens(i)/100))/sqrt(8.5))**2
	cfohe(i)=(const3*(nuo(i)+(ne(i)*hedens(i)/100))/sqrt(10.0))**2
	cfoo2(i)=(const3*(nuo(i)+(ne(i)*o2dens(i)/100))/sqrt(24.0))**2
	cfono(i)=(const3*(nuo(i)+(ne(i)*nodens(i)/100))/sqrt(23.0))**2
	cfon(i)=(const3*(nuo(i)+(ne(i)*ndens(i)/100))/sqrt(15.0))**2
	cfo2o(i)=(const3*(nuo2(i)+(ne(i)*odens(i)/100))/sqrt(24.0))**2
	cfo2h(i)=(const3*(nuo2(i)+(ne(i)*hdens(i)/100))/sqrt(16.5))**2
	cfo2he(i)=(const3*(nuo2(i)+(ne(i)*hedens(i)/100))/sqrt(18.0))**2
	cfo2o2(i)=(const3*(nuo2(i)+(ne(i)*o2dens(i)/100))/sqrt(32.0))**2
	cfo2no(i)=(const3*(nuo2(i)+(ne(i)*nodens(i)/100))/sqrt(31.0))**2
	cfo2n(i)=(const3*(nuo2(i)+(ne(i)*ndens(i)/100))/sqrt(23.0))**2
	cfn2o(i)=(const3*(nun2(i)+(ne(i)*odens(i)/100))/sqrt(22.0))**2
	cfn2h(i)=(const3*(nun2(i)+(ne(i)*hdens(i)/100))/sqrt(14.5))**2
	cfn2he(i)=(const3*(nun2(i)+(ne(i)*hedens(i)/100))/4)**2
	cfn2o2(i)=(const3*(nun2(i)+(ne(i)*o2dens(i)/100))/sqrt(30.0))**2
	cfn2no(i)=(const3*(nun2(i)+(ne(i)*nodens(i)/100))/sqrt(29.0))**2
	cfn2n(i)=(const3*(nun2(i)+(ne(i)*ndens(i)/100))/sqrt(21.0))**2
	cfheo(i)=(const3*(nuhe(i)+(ne(i)*odens(i)/100))/sqrt(21.0))**2
	cfheh(i)=(const3*(nuhe(i)+(ne(i)*hdens(i)/100))/sqrt(2.5))**2
	cfhehe(i)=(const3*(nuhe(i)+(ne(i)*hedens(i)/100))/2)**2
	cfheo2(i)=(const3*(nuhe(i)+(ne(i)*o2dens(i)/100))/sqrt(18.0))**2
	cfheno(i)=(const3*(nuhe(i)+(ne(i)*nodens(i)/100))/sqrt(17.0))**2
	cfhen(i)=(const3*(nuhe(i)+(ne(i)*ndens(i)/100))/3)**2
	cfaro(i)=(const3*(nuar(i)+(ne(i)*odens(i)/100))/sqrt(28.0))**2
	cfarh(i)=(const3*(nuar(i)+(ne(i)*hdens(i)/100))/sqrt(20.5))**2
	cfarhe(i)=(const3*(nuar(i)+(ne(i)*hedens(i)/100))/sqrt(22.0))**2
	cfaro2(i)=(const3*(nuar(i)+(ne(i)*o2dens(i)/100))/sqrt(36.0))**2
	cfarno(i)=(const3*(nuar(i)+(ne(i)*nodens(i)/100))/sqrt(35.0))**2
	cfarn(i)=(const3*(nuar(i)+(ne(i)*ndens(i)/100))/sqrt(27.0))**2
	cfho(i)=(const3*(nuh(i)+(ne(i)*odens(i)/100))/sqrt(8.5))**2
	cfhh(i)=(const3*(nuh(i)+(ne(i)*hdens(i)/100))/1.0)**2
	cfhhe(i)=(const3*(nuh(i)+(ne(i)*hedens(i)/100))/sqrt(2.5))**2
	cfho2(i)=(const3*(nuh(i)+(ne(i)*o2dens(i)/100))/sqrt(16.5))**2
	cfhno(i)=(const3*(nuh(i)+(ne(i)*nodens(i)/100))/sqrt(15.5))**2
	cfhn(i)=(const3*(nuh(i)+(ne(i)*ndens(i)/100))/sqrt(7.5))**2
	cfno(i)=(const3*(nun(i)+(ne(i)*odens(i)/100))/sqrt(15.0))**2
	cfnh(i)=(const3*(nun(i)+(ne(i)*hdens(i)/100))/sqrt(7.5))**2
	cfnhe(i)=(const3*(nun(i)+(ne(i)*hedens(i)/100))/4)**2
	cfno2(i)=(const3*(nun(i)+(ne(i)*o2dens(i)/100))/sqrt(23.0))**2
	cfnno(i)=(const3*(nun(i)+(ne(i)*nodens(i)/100))/sqrt(22.0))**2
	cfnn(i)=(const3*(nun(i)+(ne(i)*ndens(i)/100))/sqrt(14.0))**2	
	!write(*,*) nuo(i),nun2(i),nuo2(i),nuhe(i),nuar(i),nuh(i),nun(i)
	!write(1,*) height(i),cfoo(i),cfoh(i),cfohe(i),cfoo2(i),cfono(i),cfon(i),cfn2o(i),cfn2h(i),cfn2he(i),cfn2o2(i),cfn2no(i),cfn2n(i),&
	!		   cfo2o(i),cfo2h(i),cfo2he(i),cfo2o2(i),cfo2no(i),cfo2n(i),cfheo(i),cfheh(i),cfhehe(i),cfheo2(i),cfheno(i),cfhen(i),&
	!		   cfaro(i),cfarh(i),cfarhe(i),cfaro2(i),cfarno(i),cfarn(i),cfho(i),cfhh(i),cfhhe(i),cfho2(i),cfhno(i),cfhn(i),&
	!		   cfno(i),cfnh(i),cfnhe(i),cfno2(i),cfnno(i),cfnn(i)
	!! Calculation of plasma frequencies^2 follows. SI units.
	pfe(i)=c5*edens(i)/eps
	pfo(i)=const*odens(i)*edens(i)/100/16
	pfh(i)=const*hdens(i)*edens(i)/100
	pfhe(i)=const*hedens(i)*edens(i)/100/4
	pfo2(i)=const*o2dens(i)*edens(i)/100/32
	pfno(i)=const*nodens(i)*edens(i)/100/30
	pfn(i)=const*ndens(i)*edens(i)/100/14
	!write (2,*) edens(i),pfe(i),pfo(i),pfh(i),pfhe(i),pfo2(i),pfno(i),pfn(i)
	!! Calculation of gyrofrequencies^2 follows. SI units.
	gye(i)=const2*(b(i)**2)
	gyo(i)=const1*(b(i)**2)/256
	gyh(i)=const1*(b(i)**2)
	gyhe(i)=const1*(b(i)**2)/16
	gyo2(i)=const1*(b(i)**2)/1024
	gyno(i)=const1*(b(i)**2)/900
	gyn(i)=const1*(b(i)**2)/196
	!write(3,*) gye(i),gyo(i),gyh(i),gyhe(i),gyo2(i),gyno(i),gyn(i)
	!! Calculation of wave speed follows (km/s).
	alpha(i)=(pfe(i)/gye(i))+(pfo(i)*((1/(cfoo(i)+gyo(i)))+(1/(cfn2o(i)+gyo(i)))+(1/(cfo2o(i)+gyo(i)))+ &
			 (1/(cfheo(i)+gyo(i)))+(1/(cfaro(i)+gyo(i)))+(1/(cfho(i)+gyo(i)))+(1/(cfno(i)+gyo(i)))))+ &
			 (pfh(i)*((1/(cfoh(i)+gyh(i)))+(1/(cfn2h(i)+gyh(i)))+(1/(cfo2h(i)+gyh(i)))+ &
			 (1/(cfheh(i)+gyh(i)))+(1/(cfarh(i)+gyh(i)))+(1/(cfhh(i)+gyh(i)))+(1/(cfnh(i)+gyh(i)))))+ &
			 (pfhe(i)*((1/(cfohe(i)+gyhe(i)))+(1/(cfn2he(i)+gyhe(i)))+(1/(cfo2he(i)+gyhe(i)))+ &
			 (1/(cfhehe(i)+gyhe(i)))+(1/(cfarhe(i)+gyhe(i)))+(1/(cfhhe(i)+gyhe(i)))+(1/(cfnhe(i)+gyhe(i)))))+ &
			 (pfo2(i)*((1/(cfoo2(i)+gyo2(i)))+(1/(cfn2o2(i)+gyo2(i)))+(1/(cfo2o2(i)+gyo2(i)))+ &
			 (1/(cfheo2(i)+gyo2(i)))+(1/(cfaro2(i)+gyo2(i)))+(1/(cfho2(i)+gyo2(i)))+(1/(cfno2(i)+gyo2(i)))))+ &
			 (pfno(i)*((1/(cfono(i)+gyno(i)))+(1/(cfn2no(i)+gyno(i)))+(1/(cfo2no(i)+gyno(i)))+ &
			 (1/(cfheno(i)+gyno(i)))+(1/(cfarno(i)+gyno(i)))+(1/(cfhno(i)+gyno(i)))+(1/(cfnno(i)+gyno(i)))))+ &
			 (pfn(i)*((1/(cfon(i)+gyn(i)))+(1/(cfn2n(i)+gyn(i)))+(1/(cfo2n(i)+gyn(i)))+ &
			 (1/(cfhen(i)+gyn(i)))+(1/(cfarn(i)+gyn(i)))+(1/(cfhn(i)+gyn(i)))+(1/(cfnn(i)+gyn(i)))))
	speed(i)=c/((1+alpha(i))**0.5)
	va(i)=1.0E-9*b(i)/sqrt(mu*dentot(i))
	!write(*,*) i,va(i),dentot(i)
	!write(4,*) height(i),speed(i),va(i)
	!! Calculation of Pedersen conductivity follows. SI units (mho/m).
	sigped(i)=(c4/16)*(edens(i)*odens(i)/100)* &
			  ((sqrt(cfoo(i))/(cfoo(i)+gyo(i)))+(sqrt(cfn2o(i))/(cfn2o(i)+gyo(i)))+(sqrt(cfo2o(i))/(cfo2o(i)+gyo(i)))+ &
			  (sqrt(cfheo(i))/(cfheo(i)+gyo(i)))+(sqrt(cfaro(i))/(cfaro(i)+gyo(i)))+(sqrt(cfho(i))/(cfho(i)+gyo(i)))+ &
			  (sqrt(cfno(i))/(cfno(i)+gyo(i)))) + &
			  c4*(edens(i)*hdens(i)/100)* &
			  ((sqrt(cfoh(i))/(cfoh(i)+gyh(i)))+(sqrt(cfn2h(i))/(cfn2h(i)+gyh(i)))+(sqrt(cfo2h(i))/(cfo2h(i)+gyh(i)))+ &
			  (sqrt(cfheh(i))/(cfheh(i)+gyh(i)))+(sqrt(cfarh(i))/(cfarh(i)+gyh(i)))+(sqrt(cfhh(i))/(cfhh(i)+gyh(i)))+ &
			  (sqrt(cfnh(i))/(cfnh(i)+gyh(i)))) + &
			  (c4/4)*(edens(i)*hedens(i)/100)* &
			  ((sqrt(cfohe(i))/(cfohe(i)+gyhe(i)))+(sqrt(cfn2he(i))/(cfn2he(i)+gyhe(i)))+(sqrt(cfo2he(i))/(cfo2he(i)+gyhe(i)))+ &
			  (sqrt(cfhehe(i))/(cfhehe(i)+gyhe(i)))+(sqrt(cfarhe(i))/(cfarhe(i)+gyhe(i)))+(sqrt(cfhhe(i))/(cfhhe(i)+gyhe(i)))+ &
			  (sqrt(cfnhe(i))/(cfnhe(i)+gyhe(i)))) + &
			  (c4/32)*(edens(i)*o2dens(i)/100)* &
			  ((sqrt(cfoo2(i))/(cfoo2(i)+gyo2(i)))+(sqrt(cfn2o2(i))/(cfn2o2(i)+gyo2(i)))+(sqrt(cfo2o2(i))/(cfo2o2(i)+gyo2(i)))+ &
			  (sqrt(cfheo2(i))/(cfheo2(i)+gyo2(i)))+(sqrt(cfaro2(i))/(cfaro2(i)+gyo2(i)))+(sqrt(cfho2(i))/(cfho2(i)+gyo2(i)))+ &
			  (sqrt(cfno2(i))/(cfno2(i)+gyo2(i)))) + &
			  (c4/30)*(edens(i)*nodens(i)/100)* &
			  ((sqrt(cfono(i))/(cfono(i)+gyno(i)))+(sqrt(cfn2no(i))/(cfn2no(i)+gyno(i)))+(sqrt(cfo2no(i))/(cfo2no(i)+gyno(i)))+ &
			  (sqrt(cfheno(i))/(cfheno(i)+gyno(i)))+(sqrt(cfarno(i))/(cfarno(i)+gyno(i)))+(sqrt(cfhno(i))/(cfhno(i)+gyno(i)))+ &
			  (sqrt(cfnno(i))/(cfnno(i)+gyno(i)))) + &
			  (c4/14)*(edens(i)*ndens(i)/100)* &
			  ((sqrt(cfon(i))/(cfon(i)+gyn(i)))+(sqrt(cfn2n(i))/(cfn2n(i)+gyn(i)))+(sqrt(cfo2n(i))/(cfo2n(i)+gyn(i)))+ &
			  (sqrt(cfhen(i))/(cfhen(i)+gyn(i)))+(sqrt(cfarn(i))/(cfarn(i)+gyn(i)))+(sqrt(cfhn(i))/(cfhn(i)+gyn(i)))+ &
			  (sqrt(cfnn(i))/(cfnn(i)+gyn(i))))
	!! Calculation of parallel conductivity follows. SI units (mho/m)
	cfei(i)=(34+4.18*log((etemp(i)**3)/ne(i)))*ne(i)/(etemp(i)**1.5) ! electron-ion collision frequency
	sigpar(i)=c5*edens(i)/cfei(i)
	!write(5,*) height(i),sigped(i),sigpar(i)
end do

close(10)
close(11)
close(12)
close(13) 
!close(1)
!close(2)
!close(3)
!close(4)
!close(5)

!!!! Computing height-integrated Pedersen conductivity using trapezoidal rule and Simpson's rule.

hintpedz=(dz/2)*(2*sum(sigped)-sigped(1)-sigped(46))

dump1=2*sum(sigped(4:43:3))
dump2=0
dump3=0
dump4=0

do i=2,45,3
	dump2=dump2+3*sum(sigped(i:i+1))
end do

do i=4,43,3
	dump3=sigped(i)
	dump4=dump4+dump3
end do

hintpeds=(33/4)*(sigped(1)+sigped(46)+dump1+dump2)
!write(*,*) hintpeds,dump1,dump2,2*dump4
!!!!

!!!!!!!!!! Defining the coefficients multiplying the fields.

t=0
aa=dt/dz/1000 ! division by 1000 to convert to SI
ab=dt/dx/1000

!!!!!! Checking the Courant condition
!do i=1,46 
!	cour(i)=(aa*aaa(i))**0.5 
!	write(*,*) i,cour(i)
!end do
!!!!!!

! should extrapolate this code to higher altitudes by assuming an exponential or power law tail of va and all densities until 20000 km height. can expand the x direction to keep the flux through the tubes the same. should also take b out as 1/r^3.

! should take exponential(e^-r/h)+power law(1/r) for hydrogen, as it flattens out later on. 10cc/r for dentot

!!!!!!! Creating bigger arrays to go until 20000 km, i.e., z=996.

!!!!! Define an array z(i) from 1000 to 20000, 950 steps.
z(1)=1020
do i=2,950
	z(i)=1020+(i-1)*20
end do
!!!!!

!!!!! setting profiles for oxygen, hydrogen and nitrogen ion densities beyond 1000 km, and putting the rest to zero since their concentrations are already low or zero at 1000 km.
do i=1,950 
	odens1(i)=(odens(46)*edens(46)/100)*exp(-(z(i)-1000)/1000) !SI units, no longer as percentages
	hdens1(i)=hdens(46)*edens(46)*10/z(i)
	ndens1(i)=(ndens(46)*edens(46)/100)*exp(-(z(i)-1000)/1000)
	edens1(i)=odens1(i)+hdens1(i)+ndens1(i) 
	!write(6,*) i,odens1(i),hdens1(i),ndens(i),edens1(i)
end do
hedens1=0
o2dens1=0
nodens1=0

edens2(1:46)=edens
edens2(47:996)=edens1
ne1=(1.0E-6)*edens2 !converting to cm^-3
odens2(1:46)=odens*edens/100
odens2(47:996)=odens1
hdens2(1:46)=hdens*edens/100
hdens2(47:996)=hdens1
hedens2(1:46)=hedens*edens/100
hedens2(47:996)=hedens1
o2dens2(1:46)=o2dens*edens/100
o2dens2(47:996)=o2dens1
nodens2(1:46)=nodens*edens/100
nodens2(47:996)=nodens1
ndens2(1:46)=ndens*edens/100
ndens2(47:996)=ndens1

dentot1=(melec*edens2)+(mprot*16*odens2)+(mprot*hdens2)+(mprot*14*ndens2)
!!!!!

!!!!! setting all ion-neutral collision frequencies to zero beyond 1000 km as their magnitudes are negligible at that height itself.
cfoo1=0
cfoo_2(1:46)=cfoo
cfoo_2(47:996)=cfoo1

cfoh1=0
cfoh_2(1:46)=cfoh
cfoh_2(47:996)=cfoh1

cfohe1=0
cfohe_2(1:46)=cfohe
cfohe_2(47:996)=cfohe1

cfoo21=0
cfoo2_2(1:46)=cfoo2
cfoo2_2(47:996)=cfoo21

cfono1=0
cfono_2(1:46)=cfono
cfono_2(47:996)=cfono1

cfon1=0
cfon_2(1:46)=cfon
cfon_2(47:996)=cfon1

cfo2o1=0
cfo2o_2(1:46)=cfo2o
cfo2o_2(47:996)=cfo2o1

cfo2h1=0
cfo2h_2(1:46)=cfo2h
cfo2h_2(47:996)=cfo2h1

cfo2he1=0
cfo2he_2(1:46)=cfo2he
cfo2he_2(47:996)=cfo2he1

cfo2o21=0
cfo2o2_2(1:46)=cfo2o2
cfo2o2_2(47:996)=cfo2o21

cfo2no1=0
cfo2no_2(1:46)=cfo2no
cfo2no_2(47:996)=cfo2no1

cfo2n1=0
cfo2n_2(1:46)=cfo2n
cfo2n_2(47:996)=cfo2n1

cfn2o1=0
cfn2o_2(1:46)=cfn2o
cfn2o_2(47:996)=cfn2o1

cfn2h1=0
cfn2h_2(1:46)=cfn2h
cfn2h_2(47:996)=cfn2h1

cfn2he1=0
cfn2he_2(1:46)=cfn2he
cfn2he_2(47:996)=cfn2he1

cfn2o21=0
cfn2o2_2(1:46)=cfn2o2
cfn2o2_2(47:996)=cfn2o21

cfn2no1=0
cfn2no_2(1:46)=cfn2no
cfn2no_2(47:996)=cfn2no1

cfn2n1=0
cfn2n_2(1:46)=cfn2n
cfn2n_2(47:996)=cfn2n1

cfheo1=0
cfheo_2(1:46)=cfheo
cfheo_2(47:996)=cfheo1

cfheh1=0
cfheh_2(1:46)=cfheh
cfheh_2(47:996)=cfheh1

cfhehe1=0
cfhehe_2(1:46)=cfhehe
cfhehe_2(47:996)=cfhehe1

cfheo21=0
cfheo2_2(1:46)=cfheo2
cfheo2_2(47:996)=cfheo21

cfheno1=0
cfheno_2(1:46)=cfheno
cfheno_2(47:996)=cfheno1

cfhen1=0
cfhen_2(1:46)=cfhen
cfhen_2(47:996)=cfhen1

cfaro1=0
cfaro_2(1:46)=cfaro
cfaro_2(47:996)=cfaro1

cfarh1=0
cfarh_2(1:46)=cfarh
cfarh_2(47:996)=cfarh1

cfarhe1=0
cfarhe_2(1:46)=cfarhe
cfarhe_2(47:996)=cfarhe1

cfaro21=0
cfaro2_2(1:46)=cfaro2
cfaro2_2(47:996)=cfaro21

cfarno1=0
cfarno_2(1:46)=cfarno
cfarno_2(47:996)=cfarno1

cfarn1=0
cfarn_2(1:46)=cfarn
cfarn_2(47:996)=cfarn1

cfho1=0
cfho_2(1:46)=cfho
cfho_2(47:996)=cfho1

cfhh1=0
cfhh_2(1:46)=cfhh
cfhh_2(47:996)=cfhh1

cfhhe1=0
cfhhe_2(1:46)=cfhhe
cfhhe_2(47:996)=cfhhe1

cfho21=0
cfho2_2(1:46)=cfho2
cfho2_2(47:996)=cfho21

cfhno1=0
cfhno_2(1:46)=cfhno
cfhno_2(47:996)=cfhno1

cfhn1=0
cfhn_2(1:46)=cfhn
cfhn_2(47:996)=cfhn1

cfno1=0
cfno_2(1:46)=cfno
cfno_2(47:996)=cfno1

cfnh1=0
cfnh_2(1:46)=cfnh
cfnh_2(47:996)=cfnh1

cfnhe1=0
cfnhe_2(1:46)=cfnhe
cfnhe_2(47:996)=cfnhe1

cfno21=0
cfno2_2(1:46)=cfno2
cfno2_2(47:996)=cfno21

cfnno1=0
cfnno_2(1:46)=cfnno
cfnno_2(47:996)=cfnno1

cfnn1=0
cfnn_2(1:46)=cfnn
cfnn_2(47:996)=cfnn1
!!!!!

!!!!!
pfe1=c5*edens2/eps
pfo1=const*odens2/16
pfh1=const*hdens2
pfhe1=const*hedens2/4
pfo21=const*o2dens2/32
pfno1=const*nodens2/30
pfn1=const*ndens2/14
!!!!!

!!!!!
do i=1,950
	b1(i)=b(1)*((re+100)/(re+z(i)))**3
end do
b2(1:46)=b
b2(47:996)=b1
!!!!!

!!!!!
gye1=const2*(b2**2)
gyo1=const1*(b2**2)/256
gyh1=const1*(b2**2)
gyhe1=const1*(b2**2)/16
gyo21=const1*(b2**2)/1024
gyno1=const1*(b2**2)/900
gyn1=const1*(b2**2)/196
!!!!!

!!!!!
alpha1=(pfe1/gye1)+(pfo1*((1/(cfoo_2+gyo1))+(1/(cfn2o_2+gyo1))+(1/(cfo2o_2+gyo1))+ &
			 (1/(cfheo_2+gyo1))+(1/(cfaro_2+gyo1))+(1/(cfho_2+gyo1))+(1/(cfno_2+gyo1))))+ &
			 (pfh1*((1/(cfoh_2+gyh1))+(1/(cfn2h_2+gyh1))+(1/(cfo2h_2+gyh1))+ &
			 (1/(cfheh_2+gyh1))+(1/(cfarh_2+gyh1))+(1/(cfhh_2+gyh1))+(1/(cfnh_2+gyh1))))+ &
			 (pfhe1*((1/(cfohe_2+gyhe1))+(1/(cfn2he_2+gyhe1))+(1/(cfo2he_2+gyhe1))+ &
			 (1/(cfhehe_2+gyhe1))+(1/(cfarhe_2+gyhe1))+(1/(cfhhe_2+gyhe1))+(1/(cfnhe_2+gyhe1))))+ &
			 (pfo21*((1/(cfoo2_2+gyo21))+(1/(cfn2o2_2+gyo21))+(1/(cfo2o2_2+gyo21))+ &
			 (1/(cfheo2_2+gyo21))+(1/(cfaro2_2+gyo21))+(1/(cfho2_2+gyo21))+(1/(cfno2_2+gyo21))))+ &
			 (pfno1*((1/(cfono_2+gyno1))+(1/(cfn2no_2+gyno1))+(1/(cfo2no_2+gyno1))+ &
			 (1/(cfheno_2+gyno1))+(1/(cfarno_2+gyno1))+(1/(cfhno_2+gyno1))+(1/(cfnno_2+gyno1))))+ &
			 (pfn1*((1/(cfon_2+gyn1))+(1/(cfn2n_2+gyn1))+(1/(cfo2n_2+gyn1))+ &
			 (1/(cfhen_2+gyn1))+(1/(cfarn_2+gyn1))+(1/(cfhn_2+gyn1))+(1/(cfnn_2+gyn1))))

speed1=c/((1+alpha1)**0.5)
va1=1.0E-9*b2/sqrt(mu*dentot1)
!!!!!

!!!!!
!sigped_1=(c4/16)*odens2* &
!		  ((sqrt(cfoo_2)/(cfoo_2+gyo1))+(sqrt(cfn2o_2)/(cfn2o_2+gyo1))+(sqrt(cfo2o_2)/(cfo2o_2+gyo1))+ &
!		  (sqrt(cfheo_2)/(cfheo_2+gyo1))+(sqrt(cfaro_2)/(cfaro_2+gyo1))+(sqrt(cfho_2)/(cfho_2+gyo1))+ &
!		  (sqrt(cfno_2)/(cfno_2+gyo1))) + &
!		  c4*hdens2* &
!		  ((sqrt(cfoh_2)/(cfoh_2+gyh1))+(sqrt(cfn2h_2)/(cfn2h_2+gyh1))+(sqrt(cfo2h_2)/(cfo2h_2+gyh1))+ &
!		  (sqrt(cfheh_2)/(cfheh_2+gyh1))+(sqrt(cfarh_2)/(cfarh_2+gyh1))+(sqrt(cfhh_2)/(cfhh_2+gyh1))+ &
!		  (sqrt(cfnh_2)/(cfnh_2+gyh1))) + &
!		  (c4/4)*hedens2* &
!		  ((sqrt(cfohe_2)/(cfohe_2+gyhe1))+(sqrt(cfn2he_2)/(cfn2he_2+gyhe1))+(sqrt(cfo2he_2)/(cfo2he_2+gyhe1))+ &
!		  (sqrt(cfhehe_2)/(cfhehe_2+gyhe1))+(sqrt(cfarhe_2)/(cfarhe_2+gyhe1))+(sqrt(cfhhe_2)/(cfhhe_2+gyhe1))+ &
!		  (sqrt(cfnhe_2)/(cfnhe_2+gyhe1))) + &
!		  (c4/32)*o2dens2* &
!		  ((sqrt(cfoo2_2)/(cfoo2_2+gyo21))+(sqrt(cfn2o2_2)/(cfn2o2_2+gyo21))+(sqrt(cfo2o2_2)/(cfo2o2_2+gyo21))+ &
!		  (sqrt(cfheo2_2)/(cfheo2_2+gyo21))+(sqrt(cfaro2_2)/(cfaro2_2+gyo21))+(sqrt(cfho2_2)/(cfho2_2+gyo21))+ &
!		  (sqrt(cfno2_2)/(cfno2_2+gyo21))) + &
!		  (c4/30)*nodens2* &
!		  ((sqrt(cfono_2)/(cfono_2+gyno1))+(sqrt(cfn2no_2)/(cfn2no_2+gyno1))+(sqrt(cfo2no_2)/(cfo2no_2+gyno1))+ &
!		  (sqrt(cfheno_2)/(cfheno_2+gyno1))+(sqrt(cfarno_2)/(cfarno_2+gyno1))+(sqrt(cfhno_2)/(cfhno_2+gyno1))+ &
!		  (sqrt(cfnno_2)/(cfnno_2+gyno1))) + &
!		  (c4/14)*ndens2* &
!		  ((sqrt(cfon_2)/(cfon_2+gyn1))+(sqrt(cfn2n_2)/(cfn2n_2+gyn1))+(sqrt(cfo2n_2)/(cfo2n_2+gyn1))+ &
!		  (sqrt(cfhen_2)/(cfhen_2+gyn1))+(sqrt(cfarn_2)/(cfarn_2+gyn1))+(sqrt(cfhn_2)/(cfhn_2+gyn1))+ &
!		  (sqrt(cfnn_2)/(cfnn_2+gyn1)))

sigped1(1)=sigped(46)
do i=2,950
	sigped1(i)=sigped(46)*1000/z(i)
end do

sigped2(1:46)=sigped
sigped2(47:996)=sigped1
!!!!!

!!!!!
cfei1(1)=cfei(46)
do i=2,950
	cfei1(i)=cfei(46)*1000/z(i)
	!write(6,*) i,cfei1(i)
end do

cfei_2(1:46)=cfei
cfei_2(47:996)=cfei1

sigpar1=c5*edens2/cfei_2
!!!!!

!!!!!
epspar1=edens2*c5*dt*dt*10
!write (6,*) epspar1
!!!!!

!!!!!
do i=1,996
	hx(i)=sqrt(b2(1)/b2(i))
	hy(i)=sqrt(b2(1)/b2(i))
	!write (6,*) i,hx(i)
end do
!!!!!

!!!!!
rvepspar1=epspar1(996:1:-1)
beta1=alpha1(996:1:-1)
pedrev1=sigped2(996:1:-1)
rveden1=edens2(996:1:-1)
rvcfei1=cfei_2(996:1:-1)
rva1=va1(996:1:-1)
!!!!!

! tanh has a discontinuity in the first derivative, so can also use tanh^2 sometimes

!!!!!
aaa1=aa/mu/eps/(1+beta1) 
bb1=(dt/eps/(1+beta1))*pedrev1
bbb1=dt/dx/mu/rvepspar1/1000
bbb2=bbb1/hx/hy
cc1=dt*c5*rveden1
abc1=dt/rvepspar1
!!!!!

!!!!! Checking the Courant condition
cour1=(aa*aaa1)**0.5 
!!!!!

!!!!!
!do i=1,996
!	write (6,fmt='(5g12.5)') i,sigped2(i),bb1(i),bbb1(i),cour1(i)
!end do
!close(6)
!!!!!

!!!!!!! Defnining the x-grid (size 55 km) and the spatial part of the initialising functions.
x(1)=0

do i=2,12
	x(i)=x(i-1)+dx
	!write(*,*) i,x(i)
end do

do i=1,12
	ex1(i)=ex0*exp(-x(i)**2/1000)
	!write(*,*) i,ex1(i)
end do

do i=1,12 ! imposing a sheet current at the top boundary
	jdrive(i)=j0*(1-x(i)**2/1000)*exp(-x(i)**2/2000)
end do

do i=1,12 ! same as above, but with the condition imposed on By as its x-derivative is Jz
	by1(i)=j0*mu*(x(i)/1.0E-3)*exp(-x(i)**2/2000)
end do
!!!!!!!

By=0
Ex=0
Ez=0
Jz=0
Ex2=0
By2=0

do i=1,996
	Ex2(:,i)=hx(i)*Ex(:,i)
	By2(:,i)=hy(i)*By(:,i)
!	write(*,*) By2(2,i),Ex2(2,i)
end do

!!!!!
do n=1,100000
	t=t+dt
	!By2(:,1)=by1(:)*tanh(t) 
	!Jz(:,1)=jdrive(:)*tanh(t)
	Ex2(:,1)=ex1(:)*tanh(t)-rva1(1)*By2(:,1)! lets the upcoming wave fly away through the top (+va to impose another reflection)....essentially sign of the poynting vector
	!Ex2(:,996)=0 
	Ex2(:,996)=(1/mu/0.02)*By2(:,995) ! find out experimentally what values of hintped the code can tolerate. compute conductivities below 100 km and use that in hintped and see if that affects values.
	do k=1,995 ! jump by 1 now to make the k,k-1,k+1 separations work
		do i=2,12
			By2(i,k)=By2(i,k)-aa*(Ex2(i,k+1)-Ex2(i,k))+ab*(Ez(i,k)-Ez(i-1,k))
		end do
	end do ! plot Ex and Va*By and see if they're equal; they should be for a purely propagating wave with no reflections. also plot fields for various times at different z's
	do k=1,995
		do i=1,12
			Jz(i,k)=Jz(i,k)+cc1(k)*Ez(i,k)-dt*rvcfei1(k)*Jz(i,k)
		end do
	end do
	do k=1,995
		do i=1,11
			Ez(i,k)=Ez(i,k)+bbb2(k)*(By2(i+1,k)-By2(i,k))-abc1(k)*Jz(i,k)
		end do	
	end do
	do k=2,996
		do i=1,12
			Ex2(i,k)=Ex2(i,k)-aaa1(k)*(By2(i,k)-By2(i,k-1))-bb1(k)*Ex2(i,k) ! no hx in second term because Ex2 has it already
		end do
	end do
	if (mod(n,50).eq.0) then
		write(90,'(*(G0.6,:,","))') t,Ex2(2,100),By2(2,100),Jz(2,100),Ez(2,100)
		write(91,'(*(G0.6,:,","))') t,Ex2(2,200),By2(2,200),Jz(2,200),Ez(2,200)
		write(92,'(*(G0.6,:,","))') t,Ex2(2,300),By2(2,300),Jz(2,300),Ez(2,300)
		write(93,'(*(G0.6,:,","))') t,Ex2(2,400),By2(2,400),Jz(2,400),Ez(2,400)
		write(94,'(*(G0.6,:,","))') t,Ex2(2,500),By2(2,500),Jz(2,500),Ez(2,500)
		write(95,'(*(G0.6,:,","))') t,Ex2(2,600),By2(2,600),Jz(2,600),Ez(2,600)
		write(96,'(*(G0.6,:,","))') t,Ex2(2,700),By2(2,700),Jz(2,700),Ez(2,700)
		write(97,'(*(G0.6,:,","))') t,Ex2(2,800),By2(2,800),Jz(2,800),Ez(2,800)
		write(98,'(*(G0.6,:,","))') t,Ex2(2,900),By2(2,900),Jz(2,900),Ez(2,900)
		write(99,'(*(G0.6,:,","))') t,Ex2(2,980),By2(2,980),Jz(2,980),Ez(2,980)
		write(80,'(*(G0.6,:,","))') t,Ex2(5,100),By2(5,100),Jz(5,100),Ez(5,100)
		write(81,'(*(G0.6,:,","))') t,Ex2(5,200),By2(5,200),Jz(5,200),Ez(5,200)
		write(82,'(*(G0.6,:,","))') t,Ex2(5,300),By2(5,300),Jz(5,300),Ez(5,300)
		write(83,'(*(G0.6,:,","))') t,Ex2(5,400),By2(5,400),Jz(5,400),Ez(5,400)
		write(84,'(*(G0.6,:,","))') t,Ex2(5,500),By2(5,500),Jz(5,500),Ez(5,500)
		write(85,'(*(G0.6,:,","))') t,Ex2(5,600),By2(5,600),Jz(5,600),Ez(5,600)
		write(86,'(*(G0.6,:,","))') t,Ex2(5,700),By2(5,700),Jz(5,700),Ez(5,700)
		write(87,'(*(G0.6,:,","))') t,Ex2(5,800),By2(5,800),Jz(5,800),Ez(5,800)
		write(88,'(*(G0.6,:,","))') t,Ex2(5,900),By2(5,900),Jz(5,900),Ez(5,900)
		write(89,'(*(G0.6,:,","))') t,Ex2(5,980),By2(5,980),Jz(5,980),Ez(5,980)
		write(70,'(*(G0.6,:,","))') t,Ex2(8,100),By2(8,100),Jz(8,100),Ez(8,100)
		write(71,'(*(G0.6,:,","))') t,Ex2(8,200),By2(8,200),Jz(8,200),Ez(8,200)
		write(72,'(*(G0.6,:,","))') t,Ex2(8,300),By2(8,300),Jz(8,300),Ez(8,300)
		write(73,'(*(G0.6,:,","))') t,Ex2(8,400),By2(8,400),Jz(8,400),Ez(8,400)
		write(74,'(*(G0.6,:,","))') t,Ex2(8,500),By2(8,500),Jz(8,500),Ez(8,500)
		write(75,'(*(G0.6,:,","))') t,Ex2(8,600),By2(8,600),Jz(8,600),Ez(8,600)
		write(76,'(*(G0.6,:,","))') t,Ex2(8,700),By2(8,700),Jz(8,700),Ez(8,700)
		write(77,'(*(G0.6,:,","))') t,Ex2(8,800),By2(8,800),Jz(8,800),Ez(8,800)
		write(78,'(*(G0.6,:,","))') t,Ex2(8,900),By2(8,900),Jz(8,900),Ez(8,900)
		write(79,'(*(G0.6,:,","))') t,Ex2(8,980),By2(8,980),Jz(8,980),Ez(8,980)
		write(60,'(*(G0.6,:,","))') t,Ex2(11,100),By2(11,100),Jz(11,100),Ez(11,100)
		write(61,'(*(G0.6,:,","))') t,Ex2(11,200),By2(11,200),Jz(11,200),Ez(11,200)
		write(62,'(*(G0.6,:,","))') t,Ex2(11,300),By2(11,300),Jz(11,300),Ez(11,300)
		write(63,'(*(G0.6,:,","))') t,Ex2(11,400),By2(11,400),Jz(11,400),Ez(11,400)
		write(64,'(*(G0.6,:,","))') t,Ex2(11,500),By2(11,500),Jz(11,500),Ez(11,500)
		write(65,'(*(G0.6,:,","))') t,Ex2(11,600),By2(11,600),Jz(11,600),Ez(11,600)
		write(66,'(*(G0.6,:,","))') t,Ex2(11,700),By2(11,700),Jz(11,700),Ez(11,700)
		write(67,'(*(G0.6,:,","))') t,Ex2(11,800),By2(11,800),Jz(11,800),Ez(11,800)
		write(68,'(*(G0.6,:,","))') t,Ex2(11,900),By2(11,900),Jz(11,900),Ez(11,900)
		write(69,'(*(G0.6,:,","))') t,Ex2(11,980),By2(11,980),Jz(11,980),Ez(11,980)
	end if
end do
!!!!!

close(90)
close(91)
close(92)
close(93)
close(94)
close(95)
close(96)
close(97)
close(98)
close(99)
close(80)
close(81)
close(82)
close(83)
close(84)
close(85)
close(86)
close(87)
close(88)
close(89)
close(70)
close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(60)
close(61)
close(62)
close(63)
close(64)
close(65)
close(66)
close(67)
close(68)
close(69)

end program grid

