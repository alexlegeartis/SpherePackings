from fractions import Fraction as Q

C0 = Q(369, 256)
cases = (
    (10, Q(5, 128), Q(5, 4), Q(5, 4), Q(11, 40)),
    (11, Q(9, 256), Q(137, 100), Q(25, 24), Q(59, 200)),
    (12, Q(1, 27), Q(59, 40), Q(15, 16), Q(297, 1000)),
)

def lp_cap(n, q):
    s = 1 - 8*q/n
    assert s < Q(1, 2*n)
    return 2*(1-s)/(Q(1, 2*n)-s)

def sufficient(n, m, q, c):
    # If the necessary condition fails, no m-point motif has d >= q.
    pair = m*(m-1)//2
    a = 8*q/n - 1
    R = Q(m, 2) - a*pair
    if R < 0:
        return True
    B = 1-a
    beta = (C0-c)/(5*n)
    assert B >= 0 and beta > 0 and 1+a+beta*n > 0
    F0 = m + 2*(a*a*pair + (1+a)*R)
    H0 = (n*m*m-n*n*F0)/2
    test1 = H0-beta*n**3*R > max(0, 1-beta*B*n*n)*R/beta
    test2 = H0 > max(0, 1-beta*B*n*n/(1+2*beta*n))*R/beta
    return test1 or test2

for n, center, qmin, c, t0 in cases:
    factor, kappa = {10:(Q(19,20),Q(1011,1000)), 11:(Q(99,100),Q(501,500)), 12:(Q(199,200),Q(10009,10000))}[n]
    assert factor**2 * kappa**(n-1) < 1
    center *= factor
    assert 0<t0<Q(1,2)
    assert max(Q(1), n*t0*t0 + 1-2*t0) < qmin
    assert c*(qmin-(1-2*t0)) >= 1
    assert Q(113,128) <= c < C0
    top = int(lp_cap(n, qmin))
    results = {}
    for m in range(1, top+1):
        def matches(q):
            return m*m*q**n >= center**2 * 2**(2*n+2)
        if not matches(Q(2)):
            results['small'] = results.get('small',0) + 1
            continue
        if matches(qmin):
            q = qmin
        else:
            low, high = qmin, Q(2)
            for _ in range(30):
                mid = (low+high)/2
                if matches(mid):
                    high = mid
                else:
                    low = mid
            q = low
        if m > 1 and q*m*(m-1)/2 > Q(n*(m*m//4),4):
            method='cut'
        elif q > Q(n,8) and m > 2*n+1:
            method='negative'
        elif m > lp_cap(n,q):
            method='quadratic'
        elif n == 10 and q == Q(n,8) and m > 21:
            method='binary'
        else:
            assert sufficient(n, m, q, c), (n,m,q,c)
            method='frame'
        results[method] = results.get(method,0)+1
    print(n,qmin,top,dict(sorted(results.items())))

# Fixed inherited packing minimum: exact maximum motif count n <=15.
def simplex(k):
    return [tuple((u&v).bit_count()%2 for v in range(1,2**k)) for u in range(2**k)]
for n in range(2,17):
    M = (1 if n < 8 else 2 if n < 12 else 4 if n < 14 else 8 if n < 15 else 16 if n < 16 else 32)
    nextm=M+1
    if n < 16:
        assert nextm*(nextm-1) > Q(n*(nextm*nextm//4),4)
    else:
        assert sufficient(16, 33, Q(2), Q(15, 16))
    if M==1:
        words=[(0,)*n]
    elif M==2:
        words=[(0,)*n,(1,)*8+(0,)*(n-8)]
    elif M==4:
        words=[x*4+(0,)*(n-12) for x in simplex(2)]
    elif M==8:
        words=[x*2 for x in simplex(3)]
    elif M==16:
        words=simplex(4)
    else:
        words=[tuple(((u&v).bit_count()+sign)%2 for v in range(16))
               for sign in range(2) for u in range(16)]
    assert len(words)==M and all(len(x)==n for x in words)
    assert all(sum(a!=b for a,b in zip(x,y))>=8 for i,x in enumerate(words) for y in words[:i])
print('q=2 exact caps and witnesses checked for 2<=n<=16')
# Paley 12: exact witness at q=3/2.
H = [[1]*12]
for i in range(11):
    H.append([1]+[-1 if i==j else 1 if pow((j-i)%11,5,11)==1 else -1 for j in range(11)])
assert all(sum(a*b for a,b in zip(x,y))==(12 if i==j else 0)
           for i,x in enumerate(H) for j,y in enumerate(H))
W = [tuple(s*a for a in x) for s in [-1,1] for x in H]
assert len(set(W))==24
assert all(sum(a!=b for a,b in zip(x,y))>=6 for i,x in enumerate(W) for y in W[:i])
for n in range(12,17):
    assert sufficient(n,2*n+1,Q(n,8),Q(15,16))
print('Paley 12 witness and n/8 capacity controls passed')
