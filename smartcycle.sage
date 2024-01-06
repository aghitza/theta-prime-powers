# Copyright (c) 2023: Miles Koumouris and Alex Ghitza

# vi: ft=python

from sage.libs.flint.fmpz_poly import Fmpz_poly
from sage.modular.dims import sturm_bound
from sage.modular.modform.vm_basis import _delta_poly, victor_miller_basis
from sage.modular.modform.eis_series import eisenstein_series_poly, eisenstein_series_qexp
from sage.arith.all import srange
import sys

if not 'BASIS' in globals():
    BASIS = dict()

if not 'EIS' in globals():
    EIS = dict()


def _lots_of_bases_level1(kmax, verbose=False, var=None):
    prec = sturm_bound(1, kmax) + 10
    EIS[0] = Fmpz_poly(1)
    EIS[2] = Fmpz_poly(0)
    for k in [4, 6, 8, 10, 14]:
        A = eisenstein_series_poly(k, prec)
        if A[0] == -1:
            A = -A
        EIS[k] = A

    A = EIS[6]**2
    A._unsafe_mutate_truncate(prec)
    EIS[12] = A

    D = _delta_poly(prec)
    BASIS[(1, 12)] = [D] 
    for k in [0, 2, 4, 6, 8, 10, 14]:
        BASIS[(1, k)] = []
   
    if verbose:
        print("computing bases:")
    for k in srange(16, kmax, 2):
        if verbose:
            sys.stdout.write(str(k)+' ')
            sys.stdout.flush()
        prevk = k - 12
        A = EIS[prevk] * EIS[12]
        A._unsafe_mutate_truncate(prec)
        EIS[k] = A

        bas = [EIS[prevk]*D] + [f*D for f in BASIS[(1, prevk)]]
        for f in bas:
            f._unsafe_mutate_truncate(prec)
        nn = len(bas)
        for ii in range(nn):
            for jj in range(ii):
                bas[jj] = bas[jj] - bas[jj][ii+1]*bas[ii]
        BASIS[(1, k)] = bas
    if verbose:
        print()

    if var is not None:
        P = PowerSeriesRing(ZZ, var)
        if verbose:
            print("converting bases to q-expansions:")
        for k in srange(0, kmax, 2):
            if verbose:
                sys.stdout.write(str(k)+' ')
                sys.stdout.flush()
            bas = BASIS[(1, k)]
            qbas = Sequence([P(f.list()).add_bigoh(prec) for f in bas], cr=True)
            BASIS[(1, k)] = qbas
        if verbose:
            print()


def lots_of_bases(kmax, N=1, verbose=False):
    if N == 1:
        _lots_of_bases_level1(kmax, verbose)
    else:
        for w in srange(2, kmax, 2):
            if (N, w) not in BASIS:
                if verbose:
                    print("computing weight %s" %w)
                C = CuspForms(Gamma1(N), w)
                sb = C.sturm_bound() + 10
                if N == 1:
                    bas = C.q_expansion_basis(sb)
                else:
                    bas = C.q_integral_basis(sb)
                BASIS[(N, w)] = bas


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
    f = f.change_ring(R)
    if verbose:
        print(k)
    fprec = f.prec()

    # P = PowerSeriesRing(ZZ, 'q')

    while w < k:
        if verbose:
            print(w)
        sb = sturm_bound(N, w) + 10
        if sb > fprec:
            raise ValueError("need precision at least %s for the form" %sb)
        if (N, w) in BASIS:
            bas = BASIS[(N, w)]
            # if N == 1:
            #     qbas = Sequence([P(f.list()).add_bigoh(sb) for f in bas], cr=True)
            #     bas = qbas
        else:
            if N == 1:
                #bas = C.q_expansion_basis(sb)
                bas = victor_miller_basis(w, prec=sb, cusp_only=True)
            else:
                C = CuspForms(Gamma1(N), w)
                bas = C.q_integral_basis(sb)
            BASIS[(N, w)] = bas
        basis = [g.change_ring(R) for g in bas]
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
    if m > 1:
        km = 2 + 2*phipm
    else:
        km = 2 + phipm
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


def diffcycle(f, N, p, m, w, verbose=False):
    c = cycle(f, N, p, m, w, verbose)
    phipm = phi(p, m)
    if m > 1:
        tl = 2 + 2*phipm
    else:
        tl = 2 + phipm
    dc = []
    s = c[0]
    for ii in srange(len(c)-1):
        df, _ = (c[ii+1]-s-tl).quo_rem(phipm)
        dc.append(df)
        s = c[ii+1]
    df, _ = (c[0]-s-tl).quo_rem(phipm)
    dc.append(df)
    return dc


def cycper(a):
    lst = []
    for j in range(len(a)):
        lst.append(a[j:] + a[:j])
    return lst


def is_same_cycle(a, b):
    return a in cycper(b)


def eis_list(kmin, kmax, p, m, prec=3000, verbose=False):
    lst = []
    for k in srange(kmin, kmax, 2):
        if verbose:
            print(k)
        e = eisenstein_series_qexp(k, normalization="integral", prec=prec)
        dc = diffcycle(e, 1, p, m, k)
        lst.append((k, [dc.count(a) for a in range(0, min(dc)-1, -1)]))
    return lst


def cusp_onedim_list(p, m, prec=3000, verbose=False):
    lst = []
    for k in [12, 16, 18, 20, 22, 26]:
        if verbose:
            print(k)
        d = CuspForms(1, k).q_integral_basis(prec)[0]
        dc = diffcycle(d, 1, p, m, k)
        lst.append((k, [dc.count(a) for a in range(0, min(dc)-1, -1)]))
    return lst
