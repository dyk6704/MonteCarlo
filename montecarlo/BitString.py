import numpy as np             

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        return ''.join(map(str, self.config))

    def __eq__(self, other):
        if self.int() == other.int() and self.N == other.N:
            return True
        return False
    
    def __len__(self):
        return self.N 

    def on(self):
        c = 0
        for i in self.config:
            if i == 1:
                c+=1
        return c
    
    def off(self):
        c = 0
        for i in self.config:
            if i == 0:
                c+=1
        return c
    
    def flip_site(self,i):
        if self.config[i] == 0:
            self.config[i] = 1
        else:
            self.config[i] = 0
    
    def int(self):
        return int(''.join(map(str, self.config)), 2)

    def set_config(self, s):#:list[int]):
        for n in range(len(s)):
            self.config[n] = s[n]
        #self.config += s
        
    def set_int_config(self, dec:int):
        for i in range(0,self.N):
            self.config[i] = dec % 2
            dec /= 2
        self.config = np.flipud(self.config)
