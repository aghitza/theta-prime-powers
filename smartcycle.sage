# Copyright (c) 2023: Miles Koumouris and Alex Ghitza

# vi: ft=python

from sage.modular.dims import sturm_bound

BASIS = dict()

def theta(f, m=1):
    """
    Given a power series or polynomial f in one variable q, return
    q times the derivative of f with respect to q.
    """
    R = f.parent()
    q = R.gen()
    g = f
    for _ in range(m):
        g = q * g.derivative(q)
    return g


def is_in_span_of_basis(f, bas):
    """
    Return True if f is an R-linear combination of the elements of bas,
    where R is the coefficient ring of the parent of the elements of
    bas.

    Assumes that bas is echelonised.
    """
    res = sum([g * f[g.valuation()] for g in bas])
    return f == res


def phi(p,m):
    return p**(m-1)*(p-1)


def smart_filtration(f, N, p, m, k, verbose=False):
    """
    Return the p^m-filtration of a q-expansion f (p prime, m>=1) known to occur in weight k.
    We use the fact that the filtration must be congruent to k modulo phi(p^m).

    TODO: this is amenable to a binary search strategy (rather than sequentially trying all
    suitably congruent weights from the smallest possible all the way up to k).
    """
    phipm = phi(p, m)
    w = k % phipm
    if w == 0:
        w = phipm
    R = Zmod(p**m)
    if verbose:
        print(k)

    while w < k:
        if verbose:
            print(w)
        sb = sturm_bound(N, w) + 10
        if (N, p, m, w) in BASIS:
            basis = BASIS[(N, p, m, w)]
        else:
            C = CuspForms(Gamma1(N), w)
            basis = [g.change_ring(R) for g in C.q_integral_basis(sb)]
            BASIS[(N, p, m, w)] = basis
        if is_in_span_of_basis(f.truncate_powerseries(sb), basis):
            return w
        w += phipm

    return k


def cycle(f, N, p, m, w, verbose=False):
    """
    Returns the p^m-filtration theta cycle of f
    """
    R = Zmod(p**m)
    phipm = phi(p, m)
    km = 2 + 2*phipm
    k = w
    if verbose:
        print("initialising")
    g = f.change_ring(R)
    k = smart_filtration(g, N, p, m, k, verbose)
    for ii in range(m):
        if verbose:
            print(f"{ii+1} out of {m+1}")
        g = theta(g)
        k = smart_filtration(g, N, p, m, k + km, verbose)

    if verbose:
        print("cycle calculation")
    ans_list = [k]
    for ii in range(phipm-1):
        if verbose:
            print(f"{ii+1} out of {phipm}") 
        g = theta(g)
        k = smart_filtration(g, N, p, m, k + km, verbose)
        ans_list.append(k)

    return ans_list

