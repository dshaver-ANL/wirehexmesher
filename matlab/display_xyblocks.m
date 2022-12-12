%display_xyblocks

n=size(xyblocks,2);
s=0.1;
hm=size(xyblocks,5);
im=size(xyblocks,1);
jm=size(xyblocks,2);
nel=size(xyblocks,3);

size(xyblocks)


for h=1:hm
    

    clf
    axis equal
    hold on
   h 
    
    for b=1:nel
        

        %for i=1:2:im
            %for j=1:3:jm
               % plot(xyblocks(i,1:2:jm,b,1,h),xyblocks(i,1:2:jm,b,2,h),'x')
                plot(xyblocks(1:3:im,:,b,1,h),xyblocks(1:3:im,:,b,2,h),'bx')
                xx=spline(1:n,xyblocks(1,:,b,1,h),1:s:n);
                yy=spline(1:n,xyblocks(1,:,b,2,h),1:s:n);
                plot(xx,yy,'m-')
                
                hold on
                pause
            %end
           % pause
        %end
        %plot(xyblocks(im,jm,b,1,h),xyblocks(im,jm,b,2,h),'rx')
    end
    pause
end
