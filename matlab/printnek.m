% printnek.m

load nekinfo;
%load ess_xyblocks

fopen( 'data_v5f.out', 'wt' );
fopen( 'xmesh_v5f.out', 'wt' );
fopen( 'ymesh_v5f.out', 'wt' );
fopen( 'zmesh_v5f.out', 'wt' );
fopen( 'BC_v5f.out', 'wt' );
fopen( 'BCcon_v5f.out', 'wt' );

tic
dlmwrite('data_v5f.out',r_data,'precision', '%12.8f','-append');
toc
fid = fopen('xmesh_v5f.out','w')
fwrite(fid,xnek,'double')
fclose(fid)
%dlmwrite('xmesh_v5f.out',xnek,'precision', '%12.8f','-append');
toc
fid = fopen('ymesh_v5f.out','w')
fwrite(fid,ynek,'double')
fclose(fid)
%dlmwrite('ymesh_v5f.out',ynek,'precision', '%12.8f','-append');
toc
fid = fopen('zmesh_v5f.out','w')
fwrite(fid,znek,'double')
fclose(fid)
%dlmwrite('zmesh_v5f.out',znek,'precision', '%12.8f','-append');
toc
dlmwrite('BC_v5f.out', BCs, '-append');
toc
dlmwrite('BCcon_v5f.out', BCcon, '-append','precision',6 );
toc


fprintf('printnek complete \n')

%,'precision','%10.7f'
