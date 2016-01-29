function [y,ty]=myConv(x,tx,h,th)

% Esta funcion realiza la convolucion de dos señales 
% cuyas bases de tiempos son distintas de [0:length(x)-1]
% Por ejemplo, si x es una delta en n=2, podemos definirla
% como x=[1]; tx=[2]. Si h vale 2 en el intervalo [-1:1],
% se define como h=[1 1 1]; th=[-1 0 1];
%
% Uso: [y,ty]=myConv(x,tx,h,th)
%
% x,tx: señal x y su base de tiempos
% h,th: señal h y su base de tiempos
% y,ty: convolucion de x con h y su base de tiempos

nx=length(x);
nh=length(h);

fast=1;

if (fast)
    y=conv(x,h);
else
    if size(x,2)==1
       y=zeros(nx+nh-1,1);
    else
       y=zeros(1,nx+nh-1);
    end
    if nx>=nh
       % x esta fija y h se desplaza
       for i=1:nh
           y(i)=sum(x(1:i).*h(i:-1:1));
       end
       jx=2;
       hr=h(nh:-1:1);
       for i=nh+1:nx
           y(i)=sum(x(jx:i).*hr);
           jx=jx+1;
       end
       jh=2;
       jx=nx-nh+2;
       for i=nx+1:nx+nh-1
           y(i)=sum(x(jx:nx).*h(nh:-1:jh));
           jh=jh+1;
           jx=jx+1;
       end
    else
       % h esta fija y x se desplaza
       % x esta fija y h se desplaza
       for i=1:nx
           y(i)=sum(h(1:i).*x(i:-1:1));
       end
       jh=2;
       xr=x(nx:-1:1);
       for i=nx+1:nh
           y(i)=sum(h(jh:i).*xr);
           jh=jh+1;
       end
       jx=2;
       jh=nh-nx+2;
       for i=nh+1:nh+nx-1
           y(i)=sum(h(jh:nh).*x(nx:-1:jx));
           jx=jx+1;
           jh=jh+1;
       end
    end
end

tymin=tx(1)+th(1);
tymax=tx(nx)+th(nh);
if (nx>2)
    Ts=tx(2)-tx(1);
else
    Ts=th(2)-th(1); % coreccion de Alberto. Esta bien? (cambie TS por Ts)
end
ty=[0:length(y)-1]*Ts+tymin;
