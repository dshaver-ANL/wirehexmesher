function [cell2p,npc] = cell_2_pin(p2cell,ncp,ncell,npin);
%
%   Given a cell number, return a list of pins connected to it.
%

cell2p=zeros(ncell,3);npc=zeros(ncell,1);
for pin=1:npin; for j=1:ncp(pin); cell=p2cell(pin,j); 
   npc(cell)=npc(cell)+1;
   cell2p(cell,npc(cell))=pin;
end;end;
for cell=1:ncell;cell2p(cell,1:npc(cell))=sort(cell2p(cell,1:npc(cell)));end;