function [xnek,ynek,znek,BCs,BCcon]=nekify(xelems,yelems,zelems,ne)

%load 61pin_coarse5 % xyzelems, ne
% xelems(GLL R,Gll C,Gll L,Row,Column,elem,Layer)
% put all elems in one long array

Lay=size(xelems,7);
%Lay=21;
blocks=size(xelems,6);
%blocks=6;
Column=size(xelems,5);
Row=size(xelems,4);
Nx=size(xelems,1)-1;


totelem=Lay*blocks*Column*Row;
epin=(totelem/(Lay*Row*Column)) - 6*ne; % blocks  in a layer that surround pins 

% x y and z coordinates for each point
xnek=zeros(totelem,Nx+1,Nx+1,Nx+1); 
ynek=xnek;
znek=xnek;

BCs=zeros(totelem,6); % Array that gives BCs for each element
                      % 0 = no BC, 1 = Wall, 2 = Periodic
                      % sides: 1=inward, 2=outward, 3=cw, 4=ccw, 5=down, 6=up
                      % sides: 1=inward,2=cw,3=outward,4=ccw,5=down,6=up
BCcon=BCs;            % Gives corresponding element for periodic



for L=1:Lay
   for b=1:blocks
      for C=1:Column
         for R=1:Row
            elemnum=(L-1)*blocks*Column*Row+(b-1)*Column*Row+(C-1)*Row+R;
            
            % xyz coords for each elemnts GLL points
            xnek(elemnum,:,:,:)=xelems(:,:,:,R,C,b,L);
            ynek(elemnum,:,:,:)=yelems(:,:,:,R,C,b,L);
            znek(elemnum,:,:,:)=zelems(:,:,:,R,C,b,L);
            % Pin walls
            if b<(epin+1) 
               if R==1
               BCs(elemnum,4)=1;
               end  
            else         
              if R==Row
               BCs(elemnum,2)=1;
              end
            end 
            % Periodic bottom
            if L==1
               BCs(elemnum,5)=2;
               BCcon(elemnum,5)=elemnum+(Lay-1)*blocks*Column*Row;
            end
            % Periodic top
            if L==Lay
               BCs(elemnum,6)=2;
               BCcon(elemnum,6)=elemnum-(Lay-1)*blocks*Column*Row;
            end
         end
      end
   end
end


%save nekinfo xnek ynek znek BCs BCcon

% print
%{

clf
hold on
axis equal

for L=1:Lay
   for gl=1:Nx+1
      for =1:blocks
         for C=1:Column
            for R=1:Row
               for gc=1:Nx+1
                  for gr=1:Nx+1
               X=squeeze(xyzelems(1,gr,gc,gl,R,C,b,L));Y=squeeze(xyzelems(2,gr,gc,gl,R,C,b,L));
              
               %X=squeeze(xyzelems(L,b,C,R,gl,:,:,1));Y=squeeze(xyzelems(L,b,C,R,gl,:,:,2));
               %m=size(X,1);n=size(X,2); s=0.1;
               %xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
               %for i=1:m-1:m; plot(xx(i,:),yy(i,:),'m-'); end;
               %xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
               %or j=1:n-1:n; plot(xx(j,:),yy(j,:),'m-'); end; 

               plot(X,Y,'bx')
                  end
                  pause
               end
            end
         end
      end
   end
end

%}

%{
clf
hold on
axis equal

for i=1:totelem
for gl=1:Nx+1
   
      for gc=1:Nx+1
         %for gr=1:Nx+1
            X=xnek(i,:,gc,gl);
            Y=ynek(i,:,gc,gl);
            plot(X,Y,'bx')
            pause
         %end
      end
   end
end

%}

