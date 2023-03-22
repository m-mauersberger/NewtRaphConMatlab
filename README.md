# newtRaph.m

_A Matlab implementation of the Newton-Raphson method with linear constraints_


(c) Michael Mauersberger 2021 (v0.1), 2023 (v1.0), LGPL License v2.1


Newton-Raphson method with constraints for finding function roots.
Damping coefficient helps to find a feasible solution.
Number of sweep points help to set a new initial searching point on the axis between the minimum and maximum of the domain.
Search ends if function tolerance or argument tolerance has been reached.
Linear constraints are applied in the form `A*x <= b` via writing back the arguments to the constraint (done by means of QR decomposition and a pseudo-inverse).


**References:**

Wikipedia - "Newtonverfahren" (German), "Newton's method" (English)


[![View NewtRaphConMatlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://de.mathworks.com/matlabcentral/fileexchange/124745-newtraphconmatlab)


# Examples

**Damping coefficient**

```matlab
[xVal,fVal,it,ex] = newtRaph(@(x)(x.^3-2*x+2),0,[],[],[],[],1e-10,1e-10,100,.0,'cent',1e-10,2)
xVal = 8.2740e-08
fVal = 2.0000
it = 100
ex = -2

[xVal,fVal,it,ex] = newtRaph(@(x)(x.^3-2*x+2),0,[],[],[],[],1e-10,1e-10,100,.2,'cent',1e-10,2)
xVal = -1.7693
fVal = 7.6439e-11
it = 24
ex = 0
```

**Linear constraint**

```matlab
[xVal,fVal,it,ex] = newtRaph(@(x)(5*(x(1)+1)^2+x(2)^2-2*x(2)),[0 1],[],[],[],[],1e-10,1e-10,100,.2,'forw',1e-10,2)
xVal = [-0.5528  1.0000]
fVal = 8.3129e-11
it = 16
ex = 0

[xVal,fVal,it,ex] = newtRaph(@(x)(5*(x(1)+1)^2+x(2)^2-2*x(2)),[0 1],[],[],[-1 -1],[-.8],1e-8,1e-10,100,.2,'forw',1e-10,2)
xVal = [-0.5878  1.3878]
fVal = 8.1014e-09
it = 64
ex = 0
```

**Number of sweep points**

```matlab
[xVal,fVal,it,ex] = newtRaph(@(x)((x.^2-3)*(x-1)-5*exp(-100*(x+.5)^2)),0,-1,0,[],[],1e-8,1e-10,100,.2,'cent',1e-10,2)
xVal = -2.6000
fVal = 4.0000
it = 3
ex = -2

[xVal,fVal,it,ex] = newtRaph(@(x)((x.^2-3)*(x-1)-5*exp(-100*(x+.5)^2)),0,-1,0,[],[],1e-8,1e-10,100,.2,'cent',1e-10,4)
xVal = -0.4544
fVal = 7.1891e-09
it = 20
ex = 0
```
