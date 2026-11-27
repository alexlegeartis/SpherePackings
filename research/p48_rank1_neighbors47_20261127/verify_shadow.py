"""Integer q-expansion replay for the 48D shadow contradiction."""


def mul(a,b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(5)]


def power(a,k):
    z=[1,0,0,0,0]
    for _ in range(k):
        z=mul(z,a)
    return z


t3=[1,2,0,0,2]
t4=[1,-2,0,0,2]
# theta2=2*q**(1/4)*(1+q**2+q**6+...).
d8=mul([0,1,0,0,0],mul(power([1,0,1,0,0],4),power(t4,4)))
columns=[mul(power(t3,48-8*j),power(d8,j)) for j in range(4)]
a1=-columns[0][1]
a2=-columns[0][2]-a1*columns[1][2]
c4=columns[0][4]+a1*columns[1][4]+a2*columns[2][4]
assert columns[3][4] and c4%columns[3][4]==0
a3=-c4//columns[3][4]
assert (a1,a2,a3)==(-96,2400,-12700)
count=sum(c*row[3] for c,row in zip((1,a1,a2,a3),columns))
assert count==100 and count>2*48
print('PASS: forced norm-three count',count,'exceeds the orthogonal cap',2*48)
