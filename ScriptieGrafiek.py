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

#d=25309
#d=998799649

#R_d = math.log(a + d**0.5*b)
#f_1(d) = d**0.5*(math.log(4*d)+2)
#f_2(d) = math.log(2*d**0.5)

if True==False:
    import turtle
    import math

    grapher = turtle.Turtle()
    window = turtle.Screen()
    window.screensize(1100,600)
    window.title("Grafiek van de regulator")

    turtle.setup(1150, 650)

    n = 5000

    maxy = n**0.5*math.log(n)

    grapher.ht()
    grapher.speed(0)
    grapher.pu()
    grapher.goto(-500,275)

    def draw_right1(some_turtle, to_write):
        some_turtle.right(90)
        some_turtle.fd(5)

        some_turtle.pu()
        
        some_turtle.left(90)
        some_turtle.fd(5)
        some_turtle.write(to_write, align = "right")
        some_turtle.bk(5)
        some_turtle.right(90)

        some_turtle.pd()
        
        some_turtle.bk(5)
        some_turtle.left(90)

    grapher.pd()
    grapher.seth(270)

    grapher.fd(25)

    draw_right1(grapher,int(maxy+0.5))

    for i in range(4):
        grapher.fd(125)
        draw_right1(grapher,int(maxy*(3-i)/4+0.5))
        
    grapher.left(90)

    def draw_right2(some_turtle, to_write):
        some_turtle.right(90)
        some_turtle.fd(5)

        some_turtle.pu()
        
        #some_turtle.left(90)
        some_turtle.fd(15)
        some_turtle.write(to_write, align = "center")
        some_turtle.bk(15)
        #some_turtle.right(90)

        some_turtle.pd()
        
        some_turtle.bk(5)
        some_turtle.left(90)

    draw_right2(grapher,0)

    for i in range(4):
        grapher.fd(250)
        draw_right2(grapher,(i+1)*n//4)

    grapher.fd(25)


    grapher.pu()
    grapher.goto(-500,-250)
    grapher.pd()
    grapher.pencolor("red")

    for i in range(1,n+1):
        fi = i**0.5*math.log(i)
        #print(i,fi)
        grapher.goto((i/n)*1000-500,(fi/maxy)*500-250)

    grapher.pu()
    #grapher.goto(520-500,450-250)
    #grapher.write("Complexiteit bovengrens van R_d: ln(d)d^0.5",align="left", font=("Arial", 8, "normal"))


    grapher.goto(-500,-250)
    grapher.pd()
    grapher.pencolor("blue")

    for i in range(1,n+1):
        fi = math.log(2*i**0.5)
        #print(i,fi)
        grapher.goto((i/n)*1000-500,(fi/maxy)*500-250)

    grapher.pu()
    #grapher.goto(520-500,-20-250)
    #grapher.write("Complexiteit ondergrens van R_d: ln(2*d^0.5)",align="left", font=("Arial", 8, "normal"))

    grapher.pencolor("green")

    for i in range(1,n+1):
        if int(i**0.5+0.5)**2==i:
            continue
        [a,b] = smallest_solution_Pell(i)
        R_d = math.log(a+b*i**0.5)
        grapher.goto((i/n)*1000-500,(R_d/maxy)*500-250)
        grapher.dot(3,"green")

    #grapher.goto(520-500,200-250)
    #grapher.write("R_d",align="left", font=("Arial", 8, "normal"))


#################################################
#                                               #
# Numeriek bewijs voor de volgende stelling:    #
# p|y_p <=> p=2                                 #
#                                               #
#################################################

    
def primes235(limit):
    yield 2; yield 3; yield 5
    if limit < 7: return
    modPrms = [7,11,13,17,19,23,29,31]
    gaps = [4,2,4,2,4,6,2,6,4,2,4,2,4,6,2,6] # 2 loops for overflow
    ndxs = [0,0,0,0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7,7]
    lmtbf = (limit + 23) // 30 * 8 - 1 # integral number of wheels rounded up
    lmtsqrt = (int(limit ** 0.5) - 7)
    lmtsqrt = lmtsqrt // 30 * 8 + ndxs[lmtsqrt % 30] # round down on the wheel
    buf = [True] * (lmtbf + 1)
    for i in range(lmtsqrt + 1):
        if buf[i]:
            ci = i & 7; p = 30 * (i >> 3) + modPrms[ci]
            s = p * p - 7; p8 = p << 3
            for j in range(8):
                c = s // 30 * 8 + ndxs[s % 30]
                buf[c::p8] = [False] * ((lmtbf - c) // p8 + 1)
                s += p * gaps[ci]; ci += 1
    for i in range(lmtbf - 6 + (ndxs[(limit - 7) % 30])): # adjust for extras
        if buf[i]: yield (30 * (i >> 3) + modPrms[i & 7])


if True==False:
    expected = 0
    found = 0
    lim = 2
    import time
    start = time.time()
    for i in primes235(40000000):
        if i>lim:
            print(i,time.time()-start)
            lim*=2
        expected+=1/i
        [a,b] = smallest_solution_Pell(i)
        if b%i==0:
            found+=1
            print(a,b,i)

#i = 6100000
#expected =3.02
#found = 1

#i=10246417
#expected = 3.0429576393027538
#found = 1

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
        return str(self.a)+"+"+str(self.b)+"\sqrt{"+str(self.d)+"} (mod "+str(self.f)+")"

d = 5
f = 7
def func(d,f):
    ar = [[0 for j in range(f)] for i in range(f)]
    for i in range(f):
        for j in range(f):
            ar[i][j] = ZfZsqrtd(i,j,d,f)

    count = 0

    for i in range(f):
        for j in range(f):
            for i2 in range(f):
                for j2 in range(f):
                    prod = ar[i][j]*ar[i2][j2]
                    if prod.a==1 and prod.b==0:
                        count+=1
                        print(i,j)

    print(d,f,count)

def mp(a,n):
    if n==1:
        return a
    b = mp(a, n//2)
    if n% 2==0:
        return b*b
    return a*b*b

def orde(a,f):
    c = a
    count = 1
    while c.b%f != 0:
        c*=a
        count+=1
    return count
