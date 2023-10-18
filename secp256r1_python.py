import gmpy

def bi(s):
    i = 0
    for c in s:
        i <<= 8
        i |= ord(c)
    return i

def ib(i,l=32):
    s = ""
    while l:
        s = chr(0xff & i) + s
        i >>= 8
        l -= 1
    return s

# curve implementation in python
class Curve:

    def __init__(self):
        # curve parameters for NIST P-256 (ANSI prime256v1, SECG secp256r1) 
        # https://www.nsa.gov/ia/_files/nist-routines.pdf
        # http://perso.univ-rennes1.fr/sylvain.duquesne/master/standards/sec2_final.pdf
        self.p = 2**256-2**224+2**192+2**96-1
        self.a = self.p-3
        self.b = 41058363725152142129326129780047268409114441015993725554835256314039467401291
        gx = int("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16)
        gy = int("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16)
        self.g = [gx,gy]
        self.n = 115792089210356248762697446949407573529996955224135760342422259061068512044369

    def valid(self,point):
        xP = point[0]

        if xP==None:
            return False

        yP = point[1]
        return yP**2 % self.p == (pow(xP, 3, self.p) + self.a*xP + self.b) % self.p

    def decompress(self,compressed):
        byte = compressed[:2]

        # point at infinity
        if byte=="00":
            return [None,None]

        xP = int(compressed[2:], 16)
        ysqr = (pow(xP, 3, self.p) + self.a*xP + self.b) % self.p
        assert self.p % 4 == 3
        yP = pow(ysqr, (self.p + 1) // 4, self.p)
        assert pow(yP, 2, self.p)==ysqr
        if yP % 2:
            if byte=="03":
                return [xP,yP]
            return [xP, -yP % self.p]
        if byte=="02":
            return [xP,yP]
        return [xP, -yP % self.p]

    def compress(self,P):
        if P[0] == None:
            return "\x00" + "\x00"*32

        byte = "02"
        if P[1] % 2:
            byte = "03"
        return byte + hex(P[0])[2:]

    def inv(self,point):
        xP = point[0]

        if xP==None:
            return [None,None]

        yP = point[1]
        R = [xP,-yP % self.p]
        return R

    def add(self,P,Q):

        # P+P=2P
        if P==Q:
            return self.dbl(P)

        # P+0=P
        if P[0]==None:
            return Q
        if Q[0]==None:
            return P

        # P+-P=0
        if Q==self.inv(P):
            return [None,None]

        xP = P[0]
        yP = P[1]
        xQ = Q[0]
        yQ = Q[1]
        s = (yP - yQ) * gmpy.invert(xP - xQ, self.p) % self.p
        xR = (pow(s,2,self.p) - xP -xQ) % self.p
        yR = (-yP + s*(xP-xR)) % self.p
        R = [xR,yR]
        return R

    def dbl(self,P):
        # 2*0=0
        if P[0]==None:
            return P

        # yP==0
        if P[1]==0:
            return [None,None]

        xP = P[0]
        yP = P[1]
        s = (3*pow(xP,2,self.p)+self.a) * gmpy.invert(2*yP, self.p) % self.p
        xR = (pow(s,2,self.p) - 2*xP) % self.p
        yR = (-yP + s*(xP-xR)) % self.p
        R = [xR,yR]
        return R

    def mul(self, P, k):
        # x0=0
        if P[0]==None:
            return P

        N = P
        R = [None,None]

        while k:
            bit = k % 2
            k >>= 1
            if bit:
                R = self.add(R,N)
            N = self.dbl(N)

        return R

curve = Curve()
compress_pubkey = curve.compress(curve.g)
print(compress_pubkey)
print(curve.decompress(compress_pubkey))
