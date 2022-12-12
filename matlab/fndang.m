function[ang]=fndang(pin)

a=pi/3;
b=pi/2;
c=pi/6;
    
if pin==1
    add=[0 a a b a];
elseif pin==2
    add=[0 a a a b];
elseif pin==3
    add=[a a a b a];
elseif pin==4
    add=[0 a a b b];
elseif pin==7
    add=[a a a a b];
elseif pin==8
    add=[0 a b a b];
elseif pin==12
    add=[c b a a b];
elseif pin==13
    add=[0 a b b a];
elseif pin==16
    add=[c b a a a];
elseif pin==17
    add=[0 b a b a];
elseif pin==18
    add=[0 b b a a];
elseif pin==19
    add=[c a b a a];
else
    add=[0 a a a a a];    
end

if  pin>19 | pin<0
    'ohs nos'
    pin
    pause
end

ang =zeros(size(add));
ang(1)=add(1);
for i=2:size(ang,2)
    ang(i)=ang(i-1)+add(i);
end

    