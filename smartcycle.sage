# Copyright (c) 2023: Miles Koumouris

# vi: ft=python

EXP_LEN = 100

def theta(f):
    """
    Given a power series or polynomial f in one variable q, return
    q times the derivative of f with respect to q.
    """
    R = f.parent()
    q = R.gen()
    return q * f.derivative(q)


def is_in_span_of_basis(f, bas):
    """
    Return True if f is an R-linear combination of the elements of bas,
    where R is the coefficient ring of the parent of the elements of
    bas.

    Assumes that bas is echelonised.
    """

    R = bas[0].parent().base_ring()
    res = sum([g * f[g.valuation()] for g in bas])
    return f == res


def phi(p,m):
    return p**(m-1)*(p-1)

def smart_filtration(f, p, m, k):

    """
    Returns the p^m-filtration of a q-expansion f (p prime, m>=1).
    We set k to be the (nonnegative representative of) equivalence class modulo phi(p^m) of the filtration
    to speed up the calculation. If we don't know this, we set k = None.
    """
    if k==None:
        filt = 1
        jump = 1
    else:
        t = floor(k/phi(p,m))
        filt = k - t*phi(p,m)
        if filt < 2:
            filt += phi(p,m)
        jump = phi(p,m)
    in_class = False

    filt = filt - jump
    while not in_class:
        filt = filt + jump
        basis = [g.change_ring(Zmod(p ** m)) for g in CuspForms(Gamma1(N), filt).q_integral_basis(EXP_LEN)]
        if basis == []:
            continue
        in_class = is_in_span_of_basis(f.change_ring(Zmod(p ** m)), basis)

        if filt > 10000: # we shouldn't get filtrations this large for the numbers we're dealing with
            return -1
    return filt


def cycle(f, p, m, w):
    """
    Returns the p^m-filtration theta cycle of f
    """
    k = w
    ans_list = []
    i = 0
    for i in range(m):
        f = theta(f)
        k+=2
        i+=1
    while i < phi(p,m) + m:
        ans_list.append(smart_filtration(f,p,m,k))
        print(f"{round(100*(i-m) / phi(p,m),2)}%")
        f = theta(f)
        k+=2
        i+=1
    return ans_list


N = 2
w = 8
p = 7
m = 2

D = CuspForms(Gamma1(N), w).q_integral_basis(EXP_LEN)[0]
#D = ModularForms(Gamma1(N), w).q_integral_basis(EXP_LEN)[0]
print(D)




print(f"Form of weight {w} in level {N}, looking at mod {p} ^ {m} filtration theta cycle:")
# print(f"Reduction modulo {p} ^ {m}:")
# print(E)
# print(f"Theta cycle with mod {p} ^ {m} filtration:")
print(cycle(D, p, m, w))





