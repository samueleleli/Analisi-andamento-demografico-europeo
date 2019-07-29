%realizzato da Samuele Leli il 29/07/2019

clear all;
close all;
clc; 

global a b c d;

y=[547 575 601 634 656 675 692 706 721 727 728 725 735 738];

n=length(y);
npts=100;
nelem=npts*(n-1);
x=linspace(1950,2015,n);

plot(x,y,'or','MarkerSize',10);
hold on;

h=(x(n)-x(1))/nelem;

fprintf("ANDAMENTO DELLA POPOLAZIONE:\n\n");

% interpolazione tramite splines

calcola_coefficienti(n,x,y);
[yy,dyy,ddyy,xx1]=interpolazione(npts,n,x);

%trovo massimi minimi relativi e assoluti

[pmax,pmin]=trova_max_min_ass(nelem,yy);
trova_max_min_rel(npts,n,dyy,xx1,x);

%trovo l'anno in cui la popolazione ha varcato la soglia dei 700 milioni

for i=1:nelem
 if(yy(i)>700)
    soglia=xx1(i);
    break;
 end
end

%calcolo della derivata con la formula a 3 punti (aggiuntivo)
der=calcola_derivata(h,n,npts,yy);

fprintf("\nP. MAX ASS= %f \t y= %f \n",xx1(pmax),yy(pmax));
fprintf("P. MIN ASS= %f \t y= %f \n",xx1(pmin),yy(pmin));
fprintf("\nLa popolazione europea ha varcato la soglia dei 700 milioni nel %f \n",soglia);

%TASSO DI CRESCITA

fprintf("\nTASSO DI CRESCITA:\n\n");

%definisco i valori del tasso di crescita ogni 2,5 anni invece di 5
p=0;
for i=1:(2*n)-2
    p=p+npts/2;
    tasso1(i)=dyy(p)/yy(p);
end

xtasso=linspace(1950,2015,2*n-2);
%interpolazione tramite splines

calcola_coefficienti(2*n-2,xtasso,tasso1);
[tasso,dtasso,ddtasso,xx2]=interpolazione(npts,2*n-2,xtasso);

%trovo massimi,minimi relativi e assoluti

trova_max_min_rel(npts,2*n-2,dtasso,xx2,xtasso);
[pmax_tasso,pmin_tasso]=trova_max_min_ass(size(tasso,2),tasso);

fprintf("\nP. MAX ASS= %f \t y= %f \n",xx2(pmax_tasso),tasso(pmax_tasso));
fprintf("P. MIN ASS= %f \t y= %f \n",xx2(pmin_tasso),tasso(pmin_tasso));

%plot 
title('ANDAMENTO POPOLAZIONE');
plot(xx1,yy,'b');
hold on;
figure;
plot(xx1,dyy,xx1,der,'k');
hold on;
title('ANDAMENTO DERIVATA POPOLAZIONE');
figure;
plot(xtasso,tasso1,'or','MarkerSize',10);
hold on;
title('ANDAMENTO TASSO DI CRESCITA');
plot(xx2,tasso,'k');
hold on;
figure;
plot(xx2,dtasso,'b');
title('DERIVATA TASSO DI CRESCITA');
hold on;


function w=s(k,z,a,b,c,d,x)
    w=a+b*(z-x(k))+c*(z-x(k))^2+d*(z-x(k))^3;
end
function w=ds(k,z,b,c,d,x)
    w=b+2*c*(z-x(k))+3*d*(z-x(k))^2;
end
function w=dds(k,z,c,d,x)
    w=2*c+6*d*(z-x(k));
end


function calcola_coefficienti(n,x,y)
 global a b c d;
 a=y;
 dd(1)=0;
 dd(n)=0;
 A=zeros(n);
 A(1,1)=1;
 A(n,n)=1;
 for i=2:n-1
     h(i-1)=x(i)-x(i-1);
     h(i)=x(i+1)-x(i);
     f=3*(y(i+1)-y(i))/h(i)-3*(y(i)-y(i-1))/h(i-1);
     dd(i)=f;
     A(i,i)=2*(h(i-1)+h(i));
     A(i,i-1)=h(i-1);
     A(i,i+1)=h(i);
 end
 c=A\dd';
 for i=1:n-1
    d(i)=(c(i+1)-c(i))/(3*h(i));
    b(i)=(a(i+1)-a(i))/h(i)-h(i)*(2*c(i)+c(i+1))/3;
 end

end

function [yy,dyy,ddyy,xx1]=interpolazione(npts,n,x)
 global a b c d;
 cont=1;
 for k=1:n-1
    xx=linspace(x(k),x(k+1),npts);
    for i=1:npts 
        zero(cont)=0;
        xx1(cont)=xx(i);
        yy(cont)=s(k,xx(i),a(k),b(k),c(k),d(k),x);        
        dyy(cont)=ds(k,xx(i),b(k),c(k),d(k),x);
        ddyy(cont)=dds(k,xx(i),c(k),d(k),x);
        cont=cont+1;
    end
 end
end

function trova_max_min_rel(npts,n,dyy,xx1,x)
sup=1;
global a b c d ;
nmax=200;
tolass=1.0e-8;
for k=1:n-2
   inf=sup;
   sup=inf+npts;
   for i=inf:sup
     if dyy(i)*dyy(i+1)<0
      sa(1)=xx1(i);
      sb(1)=xx1(i+1); 
      for n=1:nmax
        sc(n)=0.5*(sa(n)+sb(n));
        if (abs(sb(n)-sa(n))<=tolass)
          nf=n;
          m=sc(n);
          if dyy(i)<dyy(i+1) %se derivata cresce => p.to di minimo rel
            fprintf("P. MIN REL= %f \t y= %f\tn= %d\n",m,s(k,m,a(k),b(k),c(k),d(k),x),n);
          else 
            fprintf("P. MAX REL= %f \t y= %f\tn= %d\n",m,s(k,m,a(k),b(k),c(k),d(k),x),n);
          end
          break;
        end
        fa=ds(k,sa(n),b(k),c(k),d(k),x);
        fc=ds(k,sc(n),b(k),c(k),d(k),x);
        if (fa*fc < 0.0)
            sa(n+1)=sa(n);
            sb(n+1)=sc(n);
        else
            sa(n+1)=sc(n);
            sb(n+1)=sb(n);
        end
      end
     end
   end
end
end

function der=calcola_derivata(h,n,npts,yy)

 for j=0:n-2
   p=j*npts;
   der(p+1)=(-3*yy(p+1)+4*yy(p+2)-yy(p+3))/(2*h);
    for i=(p+2):(p+npts-1)
      der(i)=(yy(i+1)-yy(i-1))/(2*h);
    end
   der(p+npts)=(yy(p+npts-2)-4*yy(p+npts-1)+3*yy(p+npts))/(2*h);
 end
end

function [pmax,pmin]=trova_max_min_ass(nelem,yy)
    pmax=1;
    pmin=1;
    for i=1:nelem
        if yy(i)>yy(pmax)
            pmax=i;
        end
        if yy(i)<yy(pmin)
            pmin=i;
        end
    end
end