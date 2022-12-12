% displayelems.m
fprintf('Running displayelems.m \n')


% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whole_ind = 0; % 0 to show whole layer, 1 to show individual pins

points  = 0; % Toggle to display xs on gll points
splines = 1; % Toggle to display lines connecting gll points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load 61pin_v5% pin_blcko, pin_block, type_block, xelems, yelems, zelems, ne
% xelems(GLL R, GLL C, GLL L, Row, Column, elem, Layer)


Lay=size(xelems,7);
elem=size(xelems,6);
Col=size(xelems,5);
Row=size(xelems,4);
Nx=size(xelems,1)-1;


if whole_ind==1
    elem=size(pin_blcko,2)-ne*6;
    e_use=zeros(elem,1);
    for i=1:elem
        e_use(i)=i+(pin_blcko(i)-pin_block(i))*6;
    end
else
    e_use=1:elem;
end


clf
hold on
axis equal

%{
Lay=1;
 for gl=1:Nx+1
    e=114;
        for C=1:Col
            for R=1:Row
               X=squeeze(xelems(:,:,gl,R,C,e,L));
               Y=squeeze(yelems(:,:,gl,R,C,e,L));
               plot(X,Y,'bx')
            end
        end
        pause
        clf
        hold on
        axis equal
 end

disp "here"
pause
%}
ec = 0;
pin_old=1;
for L=15:Lay
    L

    figure(1)
    clf
    hold on
   for gl=1:2:Nx

      for e=1:elem 
      %for e=338:342    
         for C=1:Col
            for R=1:Row
                
             ec = ec+1

%               if (ec>=769-100)
 
               X=squeeze(xelems(:,:,gl,R,C,e_use(e),L));
               Y=squeeze(yelems(:,:,gl,R,C,e_use(e),L));
               %X=squeeze(xelems(:,:,gl,R,C,e,L));
               %Y=squeeze(yelems(:,:,gl,R,C,e,L));

               if splines==1
                  m=size(X,1);n=size(X,2); s=0.2;
                  xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
                  for i=1:m-1:m; plot(xx(i,:),yy(i,:),'m-'); end;
                  xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
                  for j=1:n-1:n; plot(xx(j,:),yy(j,:),'m-'); end; 
               end


               if points==1
                  plot(X,Y,'b.')
               end
              
               if (ec>5300)                
               pause    
               end
%               end  % endif ec               


            end
         end

         if whole_ind==1 && pin_blcko(e)~=pin_blcko(e+1)
                                       set(gca, 'visible', 'off') ;
            set(gcf,'color','w');
           axis equal 
             pause
             clf

             hold on
             axis equal

            
         end
         
      end
      
      
      %%%%%%%%%color
      %{
      for e=71:76%elem 
         for C=1:Col
            for R=1:Row
                
                
                
               X=squeeze(xelems(:,:,gl,R,C,e_use(e),L));
               Y=squeeze(yelems(:,:,gl,R,C,e_use(e),L));
               %X=squeeze(xelems(:,:,gl,R,C,e,L));
               %Y=squeeze(yelems(:,:,gl,R,C,e,L));
               
               if splines==1
                  m=size(X,1);n=size(X,2); s=0.2;
                  xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
                  for i=1:m-1:m; plot(xx(i,:),yy(i,:),'b-'); end;
                  xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
                  for j=1:n-1:n; plot(xx(j,:),yy(j,:),'b-'); end; 
               end


               if points==1
                  plot(X,Y,'bx')
               end
        
                
            end
         end
         %pause
         if whole_ind==1 && pin_blcko(e)~=pin_blcko(e+1)
                                       set(gca, 'visible', 'off') ;
            set(gcf,'color','w');
             pause
             clf

             hold on
             axis equal

             
         end
      end 
      %end color
            %%%%%%%%%color
      for e=71%elem 
         for C=1:Col
            for R=1:Row
                
                
                
               X=squeeze(xelems(:,:,gl,R,C,e_use(e),L));
               Y=squeeze(yelems(:,:,gl,R,C,e_use(e),L));
               %X=squeeze(xelems(:,:,gl,R,C,e,L));
               %Y=squeeze(yelems(:,:,gl,R,C,e,L));
               
               if splines==1
                  m=size(X,1);n=size(X,2); s=0.2;
                  xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
                  for i=1:m-1:m; plot(xx(i,:),yy(i,:),'m-'); end;
                  xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
                  for j=1:n-1:n; plot(xx(j,:),yy(j,:),'m-'); end; 
               end


               if points==1
                  plot(X,Y,'bx')
               end
        
                
            end
         end
         %pause
         if whole_ind==1 && pin_blcko(e)~=pin_blcko(e+1)
                                       set(gca, 'visible', 'off') ;
            set(gcf,'color','w');
             pause
             clf

             hold on
             axis equal

             
         end
      end 
      %end color
      %}
      
      if whole_ind==0
            set(gca, 'visible', 'off') ;
            set(gcf,'color','w');
         pause
         clf
         hold on
         axis equal
      end
   end
end
fprintf('finished')


