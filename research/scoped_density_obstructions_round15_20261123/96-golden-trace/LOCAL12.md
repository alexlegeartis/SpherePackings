# Local contact cap in the 96D norm-12 golden shell

In the cancellation setting of [README.md](README.md), let \(X\) be the
norm-12 shell with \(s=0\). For any \(x\in X\), fewer than 313,000,000
vectors \(y\in X\) have \(\langle x,y\rangle=8\). This is a bound on the
nearest layer inside one shell, not on the whole shell.

Map such \(y\) injectively to \(r=x-y\). Then
\[
 r^2=8,\quad \langle r,Sr\rangle=0,\quad
 \langle r,x\rangle=4,\quad \langle r,Sx\rangle=0.
\]
For the last identity, \(s(x-y)=-2\langle x,Sy\rangle=0\) by
shortest-shell cancellation. The vectors \(x\) and \(Sx\) are nonzero
and orthogonal since \(s(x)=0\). Subtracting \(x/3\) from every \(r\)
leaves a vector in their common orthogonal complement, of dimension 94
and squared norm \(20/3\). Distinct \(r,r'\) cannot be opposites since
both have dot product 4 with \(x\). The parent minimum then gives
\(\langle r,r'\rangle=j\in\{-4,\ldots,4\}\), since both
\(r+r'\) and \(r-r'\) have squared norm at least eight and the parent
is integral. The resulting unit vectors therefore have off-diagonal
inner products \((3j-4)/20\).

The recurrence and eight positive rational coefficients in local12.py
give a Delsarte polynomial with constant normalized Gegenbauer
coefficient 1 in dimension 94. All nine off-diagonal values are negative
and its diagonal value is less than 313,000,000. The addition formula
bounds the number of unit vectors by that diagonal value. The two
signs \(\langle x,y\rangle=\pm8\) have equal counts by antipodality;
their combined count is also at most twice the norm-8 shell size.
This gives a multiplicity constraint for further shell estimates.

Run python -B local12.py for an exact standard-library replay.
