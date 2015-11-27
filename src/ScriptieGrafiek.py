import math
import matplotlib.pyplot as plt

def isprime(n):
    for i in range(2,math.ceil(math.sqrt(n))+1):
        if n%i==0:
            return False
    return True

def modinv(n,m):
    """Finds an x with nx + my = 1 (the inverse of n mod m), assuming that ggd(n,m)=1"""
    x_prev = 1
    x = 0
    ggd_prev = n
    ggd = abs(m)
    while ggd!=1:
        ratio = ggd_prev//ggd
        x_prev,x = x,x_prev-ratio*x
        ggd_prev,ggd = ggd,ggd_prev%ggd
    return x%m

def smallest_solution_Pell(d):
    """Finds the smallest positive integers a,b that fulfill Pell's equation: x**2-d*y**2=1, using the Chakravala method"""
    a,b = int(d**0.5)+1,1
    k = a**2-d*b**2
    #print(a,b,k)
    while k!=1:
        m_0 = (-a*modinv(b,k))%k
        m_1 = (int(d**0.5+0.5)//k)*k + m_0
        p_1 = m_1**2-d
        if p_1>0:
            m_2 = m_1 - abs(k)
        else:
            m_2 = m_1 + abs(k)
        p_2 = m_2**2-d
        if abs(p_1)<abs(p_2):
            m = m_1
        else:
            m = m_2
        #print(a,b,k,m_0,m_1,p_1,m_2,p_2,m)
        a,b,k = (a*m+d*b)//abs(k),(a+b*m)//abs(k),(m**2-d)//k
        #print(a,b,k,a**2-d*b**2)
    return [a,b]

class ZfZsqrtd():
    def __init__(self,a,b,d,f=-1):
        if f==-1:
            self.a = a
            self.b = b
        else:
            self.a = a%f
            self.b = b%f
        self.d = d
        self.f = f
    def __add__(self,other):
        if self.d == other.d and self.f == other.f:
            return ZfZsqrtd(self.a + other.a, self.b + other.b, self.d, self.f)
    def __mul__(self, other):
        if self.d == other.d and self.f == other.f:
            return ZfZsqrtd(self.a*other.a + self.d*self.b*other.b, self.a*other.b + self.b*other.a, self.d, self.f)
    def __str__(self):
        if f!=-1:
            return str(self.a)+"+"+str(self.b)+"sqrt{"+str(self.d)+"} (mod "+str(self.f)+")"
        return str(self.a)+"+"+str(self.b)+"sqrt{"+str(self.d)+"}"

def orde(a,f):
    c = a
    count = 1
    while c.b%f != 0:
        c*=a
        count+=1
    return count

def graph(points,Xscale=1000,logaritmic=False):
    fig=plt.figure()
    plot=fig.add_subplot(111)
    x=[]
    y=[]
    prime=[]

    for i in range(2,points):
        try:
            sol=smallest_solution_Pell(i)
            x.append(sol[0])
            y.append(sol[1])
        except:
            x.append(0)
            y.append(0)    
        if isprime(i):
            prime.append(1)
        else:
            prime.append(0)
    for (xc,yc,n,p) in zip(x,y,range(2,points), prime):
        try:
            plot.scatter(xc,yc,s=10,color=(n/points,0,p))
        except:
            'do nothing'
    if logaritmic:
        plot.set_yscale('log')
        plot.set_xscale('log')
    plt.axis([1, 10*Xscale, 1, Xscale])
    plt.show()   
