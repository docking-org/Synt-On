import sys, os
from xmlrpc.client import ServerProxy
srcPath = os.path.split(os.path.realpath(__file__))[0]

if __name__ == '__main__':     
    port = sys.argv[1]
    funct = sys.argv[2]
    lib = sys.argv[3]
    out = sys.argv[4]
    
    proxy = ServerProxy("http://localhost:" + port, allow_none=True)
        
    if funct == 'analog':
        proxy.analog(lib, os.path.join(srcPath, out))