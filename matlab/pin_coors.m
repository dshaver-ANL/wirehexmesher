function [npin,X,Y,xm,ym,km] = pin_coors(ne,P);

%% ne = # pins on an edge  (2,3,...,etc.)
%% P  = pin pitch


nr = 2*ne-1;                   % # pin rows and max in any row
xm=zeros(nr,nr); ym=xm; km=xm; % Matrix of pin coordinates

dxr=P/2; dyr=dxr*sqrt(3.);
x0=0; y0=0; k=0;
nc = ne;
for ir=1:nr;
    for jc=1:nc;
       x=x0+(jc-1)*P; y=y0+dyr*(ir-1); k=k+1;
%      txtcirc(x,y,D,sprintf('%d',k)); axis equal; hold on;
       X(k)=x; Y(k)=y; xm(ir,jc)=x; ym(ir,jc)=y; km(ir,jc)=k;
    end;
    if ir<ne; nc=nc+1; x0=x0-dxr; elseif ir<nr; nc=nc-1; x0=x0+dxr; end;
end;

npin=k; x0=sum(X)/npin; y0=sum(Y)/npin; X=X-x0;Y=Y-y0;xm=xm-x0;ym=ym-y0;

