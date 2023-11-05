# Compute prime-power theta cycles of modular forms
SageMath code accompanying Miles Koumouris's MSc thesis.

Run as follows:
```
sage: %attach smartcycle.sage
sage: f = CuspForms(Gamma1(2), 8).q_integral_basis(1000)[0]
sage: %time c = cycle(f, 2, 5, 2, 8, False)
```
After a little while you get
```
CPU times: user 17.8 s, sys: 47.9 ms, total: 17.9 s
Wall time: 18 s
sage: c
[52, 34, 56, 38, 60, 62, 44, 46, 48, 50, 32, 54, 56, 58, 40, 42, 44, 66, 48, 50]
```
