import numpy as np
import math
import os
import time

def main():
    x = np.fromfile('xmesh_v5f.out',dtype=float)
    x = np.reshape(x,(-1,27),order='F')
    np.savetxt('xmesh_v5f.out',x,fmt='%12.8f',delimiter=',')

    print("printed in %s"%(time.time()-starttime))
    
    y = np.fromfile('ymesh_v5f.out',dtype=float)
    y = np.reshape(y,(-1,27),order='F')
    np.savetxt('ymesh_v5f.out',y,fmt='%12.8f',delimiter=',')

    print("printed in %s"%(time.time()-starttime))
    
    z = np.fromfile('zmesh_v5f.out',dtype=float)
    z = np.reshape(z,(-1,27),order='F')
    np.savetxt('zmesh_v5f.out',z,fmt='%12.8f',delimiter=',')

    print("printed in %s"%(time.time()-starttime))
    
    return


if __name__=="__main__":
    starttime = time.time()
    main()
    print('--- Code ran in %s seconds ---'%(time.time()-starttime))
