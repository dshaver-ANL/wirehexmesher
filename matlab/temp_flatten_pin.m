

% Playing with geometry of wires to flatten near approach of pin

rp1=6.55/2;
rp2=6.55/2;
rw1=1.85/2;

cp1=[0 0];
cp2=[8.4 0];
cwr=(6.55+1.85)/2;
limit=5.025;



t=-1:0.002:1;



pin1=[rp1*cos(pi*t)+cp1(1); rp1*sin(pi*t)+cp1(2)]';
pin2=[rp2*cos(pi*t)+cp2(1); rp2*sin(pi*t)+cp2(2)]';
cw1=[cwr*cos(pi*t)+cp1(1); cwr*sin(pi*t)+cp1(2)]';




for i=1:size(t')-1

  
    clf   
    wire1=[rw1*cos(pi*t)+cw1(i,1); rw1*sin(pi*t)+cw1(i,2)]';
    
    if any(wire1(:,1)>limit)
        for j=1:size(wire1)
            if wire1(j,1)>limit
                wire1(j,1)=limit;
            end
        end
    end
    
    
    hold on
    plot(pin1(:,1),pin1(:,2))
    plot(pin2(:,1),pin2(:,2))
    plot(wire1(:,1),wire1(:,2))
    axis equal
    
    pause
end


