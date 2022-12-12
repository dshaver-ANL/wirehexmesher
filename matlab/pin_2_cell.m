%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Couple cells to pins and generate pin-based mesh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p2cell,ncp] = pin_2_cell(nvtx,ntcl,tcell,necl,ecell,nccl,ccell);
p2cell = zeros(nvtx,7); ncp=zeros(nvtx,1);

for j=1:ntcl; for i=1:3; pin=tcell(i,j); 
   ncp(pin)=ncp(pin)+1; p2cell(pin,ncp(pin))=j; end;end;

for j=1:necl; for i=1:4; pin=ecell(i,j); k=j+ntcl;
   ncp(pin)=ncp(pin)+1; p2cell(pin,ncp(pin))=k; end;end;

for j=1:nccl; for i=1:4; pin=ccell(i,j); k=j+ntcl+necl;
   ncp(pin)=ncp(pin)+1; p2cell(pin,ncp(pin))=k; end;end;

