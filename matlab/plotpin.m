nnn = 202390
npin = 184320

lz1 = 1

figure(300)
clf
for ii = 1:elemnum
if (abs(ii- (nnn-npin) )<20)
ii
figure(300)
hold on
plot(squeeze(pmeshx(ii,1,:,lz1)),squeeze(pmeshy(ii,1,:,lz1)),'mx')
plot(squeeze(pmeshx(ii,2:3,:,lz1)),squeeze(pmeshy(ii,2:3,:,lz1)),'kx')

if (ii==(nnn-npin))
plot(squeeze(pmeshx(ii,1,:,lz1)),squeeze(pmeshy(ii,1,:,lz1)),'gx')
plot(squeeze(pmeshx(ii,2:3,:,lz1)),squeeze(pmeshy(ii,2:3,:,lz1)),'gx')
end





pause
end
end
%}


