  function [Dprim, Soper, awzp] = Helsing
% *** program for circle and star (sign correction/secant method) ***
%   close all
%   clear
  format long
%
  T16=Tinit16;
  W16=Winit16;
  npan=100;
  Radius=2;
  [z,zp,zpp,nz,~,zinter,~,w,awzp,np]=zinit(Radius,npan,T16,W16);
  Dprim=-M1Ainit(z,zp,zpp,nz,w,awzp)/pi;
  Soper=MLinitC(z,nz,awzp,zinter,npan,1)/pi;
%
% *** error in eigenvalue number "tpt" ***
%  tpt=210;
%  lambda_correct=correct(tpt)
%  lambda_fast=lambda(tpt)
%  lambda_fast_error=abs(lambda_fast-lambda_correct)
%
% ***  ***
  
  function myeig=geteig(Dprim,Soper,lambdaguess,np)
  myeig=real(eigs(eye(np)-Dprim+lambdaguess*Soper,1,'sm'));
  
  function [z,zp,zpp,nz,s,zinter,wzp,w,awzp,np]=zinit(R,npan,T16,W16)
  sinter=linspace(-pi,pi,npan+1);
  sinterdiff=ones(npan,1)*(2*pi/npan);
  a=0;
  ns=0;
  np=16*npan
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*16+1:k*16;
    sdif=sinterdiff(k)/2;
    s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T16;
    w(myind)=W16*sdif;
  end
  z=zfunc(s,R,a,ns);
  zp=zpfunc(s,R,a,ns);
  zpp=zppfunc(s,R,a,ns);
  zinter=zfunc(sinter,R,a,ns);
  nz=-1i*zp./abs(zp);
  wzp=w.*zp;
  awzp=abs(wzp);
  
  function zout=zfunc(p,R,a,ns)
  zout=R*(1+a*cos(ns*p)).*exp(1i*p);

  function zpout=zpfunc(p,R,a,ns)
  zpout=R*(1i-ns*a*sin(ns*p)+1i*a*cos(ns*p)).*exp(1i*p);   

  function zppout=zppfunc(p,R,a,ns)   
  zppout=-R*(1+(1+ns*ns)*a*cos(ns*p)+2*1i*ns*a*sin(ns*p)).*exp(1i*p);   

  function M1=M1Ainit(z,zp,zpp,nz,w,awzp)
% *** adjoint of double layer potential ***   
  N=length(z); 
  M1=zeros(N);
  for n=1:N
    M1(:,n)=awzp(n)*real(nz./(z(n)-z));
  end
  M1(1:N+1:N^2)=-w.*imag(zpp./zp)/2;
  
  function M1=MLinitC(z,nz,awzp,zinter,nseg,iper)
% *** Logarithmic potential log(|tau-z|)dt ***   
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
% *** samma metod som när målpunkterna ligger utanför randen ***
  N=length(z);
  Ng=16;
  M1=zeros(N);
  for m=1:N
    M1(:,m)=log(abs(z(m)-z))*awzp(m);
  end
  if iper==1
    kstart=1;
    kend=nseg;  
  else
    kstart=2;
    kend=nseg-1;  
  end
% *** central blocks ***
  for k=1:nseg
    myind=k*16-15:k*16;
    TEST=LGIcompRecFASon(z(myind),z(myind),zinter(k),zinter(k+1), ...
			 nz(myind),awzp(myind));
    M1(myind,myind)=TEST.*awzp(myind,ones(1,Ng)).';
  end
% *** superdiagonal blocks (targets to the left) ***
  for k=kstart:nseg
    myinds=(k-1)*Ng+1:k*Ng;
    km1=mod(k-2,nseg)+1;
    myindt=(km1-1)*Ng+1:km1*Ng;       
    [TEST,accept]=LGIcompRecFASon(z(myindt),z(myinds),zinter(k), ...
				  zinter(k+1),nz(myinds),awzp(myinds));
    na=length(accept);
    M1(myindt(accept),myinds)=TEST.*awzp(myinds,ones(1,na)).';
  end
% *** subdiagonal blocks (targets to the right) ***
  for k=1:kend
    myinds=(k-1)*Ng+1:k*Ng;
    kp1=mod(k,nseg)+1;
    myindt=(kp1-1)*Ng+1:kp1*Ng;       
    [TEST,accept]=LGIcompRecFASon(z(myindt),z(myinds),zinter(k), ...
				  zinter(k+1),nz(myinds),awzp(myinds));
    na=length(accept);
    M1(myindt(accept),myinds)=TEST.*awzp(myinds,ones(1,na)).';
  end
  
  function [LGIV,accept]=LGIcompRecFASon(ztarg,zsource,za,zb,nz,awzp)
% *** absolute values, transformed ztarg and z ***
  zmid=zb+za;
  zdif=zb-za;     
  cc=zdif/2;
  ztgtr=(2*ztarg-zmid)/zdif;
  zsctr=(2*zsource-zmid)/zdif;
  Ng=length(zsctr);
  tol=2.0;
  accept=1:Ng;
  accept=accept(abs(ztgtr)<tol);
  alen=length(accept);
  A=ones(Ng);
  for k=2:Ng
    A(:,k)=A(:,k-1).*zsctr;
  end
  p=zeros(Ng+1,1);
  Q=zeros(Ng,alen);
  c=((1-(-1).^(1:Ng))./(1:Ng))';
  for jj=1:alen
    j=accept(jj);
    upp=log(1-ztgtr(j));
    loo=log(-1-ztgtr(j));
    p(1)=upp-loo;
    p111=upp+loo;
    for k=2:Ng+1
      p(k)=ztgtr(j)*p(k-1)+c(k-1);
    end
    Q(1:2:Ng-1,jj)=(p111-p(2:2:Ng))./(1:2:Ng-1)';
    Q(2:2:Ng,jj)=(p(1)-p(3:2:Ng+1))./(2:2:Ng)';
  end
  tmp=cc*conj(nz)./awzp;
  LGIV=imag(1i*log(cc)+Q.'/A.*tmp(:,ones(1,alen)).');

  function T=Tinit16
% *** 16-point Gauss-Legendre nodes ***  
  T=zeros(16,1);
  T( 1)=-0.989400934991649932596154173450332627;
  T( 2)=-0.944575023073232576077988415534608345;
  T( 3)=-0.865631202387831743880467897712393132;
  T( 4)=-0.755404408355003033895101194847442268;
  T( 5)=-0.617876244402643748446671764048791019;
  T( 6)=-0.458016777657227386342419442983577574;
  T( 7)=-0.281603550779258913230460501460496106;
  T( 8)=-0.095012509837637440185319335424958063;
  T( 9)= 0.095012509837637440185319335424958063;
  T(10)= 0.281603550779258913230460501460496106;
  T(11)= 0.458016777657227386342419442983577574;
  T(12)= 0.617876244402643748446671764048791019;
  T(13)= 0.755404408355003033895101194847442268;
  T(14)= 0.865631202387831743880467897712393132;
  T(15)= 0.944575023073232576077988415534608345;
  T(16)= 0.989400934991649932596154173450332627;

  function W=Winit16
% *** 16-point Gauss-Legendre weights ***  
  W=zeros(16,1); 
  W( 1)= 0.027152459411754094851780572456018104;
  W( 2)= 0.062253523938647892862843836994377694;
  W( 3)= 0.095158511682492784809925107602246226;
  W( 4)= 0.124628971255533872052476282192016420;
  W( 5)= 0.149595988816576732081501730547478549;
  W( 6)= 0.169156519395002538189312079030359962;
  W( 7)= 0.182603415044923588866763667969219939;
  W( 8)= 0.189450610455068496285396723208283105;
  W( 9)= 0.189450610455068496285396723208283105;
  W(10)= 0.182603415044923588866763667969219939;
  W(11)= 0.169156519395002538189312079030359962;
  W(12)= 0.149595988816576732081501730547478549;
  W(13)= 0.124628971255533872052476282192016420;
  W(14)= 0.095158511682492784809925107602246226;
  W(15)= 0.062253523938647892862843836994377694;
  W(16)= 0.027152459411754094851780572456018104;
  