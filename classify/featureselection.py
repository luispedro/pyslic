from __future__ import division
from numpy import *
from numpy.linalg import det
from classifier import normaliselabels
import scipy.stats
import warnings

TOLERANCE = 0
SIGNIFICANCE_IN = .15
SIGNIFICANCE_OUT = .15

def _sweep(A,k,flag):
    N,_=A.shape
    Akk=A[k,k]
    B=zeros_like(A)
    try:
        from scipy import weave
        from scipy.weave import converters
        k=int(k)
        code = '''
#line 32 "featureselection.py"
        for (int i = 0; i != N; ++i) {
            for (int j = 0; j != N; ++j) {
                if (i == k) {
                    if (j == k) {
                        B(i,j) =  - 1./A(k,k);
                    } else {
                        B(i,j)=flag*A(i,j)/A(k,k);
                    }
                } else if (j == k) {
                    B(i,j)=flag*A(i,j)/A(k,k);
                } else { 
                    B(i,j)=A(i,j) - A(i,k)*A(k,j)/A(k,k);
                }
            }
        }
        '''
        weave.inline(
                code,
                ['A','B','k','Akk','flag','N'],
                type_converters=converters.blitz)
    except:
        for i in xrange(N):
            for j in xrange(N):
                if i == k:
                    if j == k:
                        B[i,j] =  - 1./Akk
                    else:
                        B[i,j]=flag*A[i,j]/Akk
                elif j == k:
                    B[i,j]=flag*A[i,j]/Akk
                else:
                    B[i,j]=A[i,j] - A[i,k]*A[k,j]/Akk
    return B

def sda(features,labels):
    '''
    features_idx = sda(features,labels)

    Perform Stepwise Discriminant Analysis for feature selection
    '''


    N, m = features.shape
    labels,labelsu = normaliselabels(labels)
    q=len(labelsu)

    # This is how the code in ml_stepdisc computes F_in (F_out)

    mus=array([features[labels==i,:].mean(0) for i in xrange(q)])
    mu=features.mean(0)
    
    W=zeros((m,m))
    T=zeros((m,m))
    try:
        from scipy import weave
        from scipy.weave import converters
        code='''
#line 68 "featureselection.py"
        for (int i = 0; i != m; ++i) {
            for (int j = 0; j != m; ++j) {
                for (int n = 0; n != N; ++n) {
                    int g=labels(n);
                    W(i,j) += (features(n,i)-mus(g,i))*(features(n,j)-mus(g,j));
                    T(i,j) += (features(n,i)-mu(i))*(features(n,j)-mu(j));
                }
            }
        }
        '''
        weave.inline(
                code,
                ['N','m','W','T','features','mu','mus','labels'],
                type_converters=converters.blitz)
    except:
        warnings.warn('scipy.weave failed. Resorting to (slow) Python code')
        for i in xrange(m):
            for j in xrange(m):
                for n in xrange(N):
                    g=labels[n]
                    W[i,j] += (features[n,i]-mus[g,i])*(features[n,j]-mus[g,j])
                    T[i,j] += (features[n,i]-mu[i])*(features[n,j]-mu[j])
    ignoreidx = ( W.diagonal() == 0 )
    if ignoreidx.any():
        F=sda(features[:,~ignoreidx],labels)
        return F + cumsum(ignoreidx)

    output=[]
    D=W.diagonal()
    df1 = q-1
    while True:
        V = W.diagonal()/T.diagonal() 
        W_d = W.diagonal()
        V_neg = (W_d < 0)
        p=V_neg.sum()
        df2 = N-p-q+1
        if V_neg.any(): 
            V_m = V[V_neg].min()
            k,=where(V == V_m)
            k=k[0]
            Fremove = (N-p-q+1)/(q-1)*(V_m-1)
            PrF = 1 - scipy.stats.f.cdf(Fremove,df1,df2)
            if PrF > SIGNIFICANCE_OUT:
                #print 'removing ',k, 'V(k)', 1./V_m, 'Fremove', Fremove
                W=_sweep(W,k,1)
                T=_sweep(T,k,1)
                continue
        ks = ( (W_d / D) > TOLERANCE)
        if ks.any():
            V_m=V[ks].min()
            k,=where(V==V_m)
            k=k[0]
            Fenter = (N-p-q)/(q-1) * (1-V_m)/V_m
            PrF = 1 - scipy.stats.f.cdf(Fenter,df1,df2)
            if PrF < SIGNIFICANCE_IN:
                #print 'adding ',k, 'V(k)', 1./V_m, 'Fenter', Fenter
                W=_sweep(W,k,-1)
                T=_sweep(T,k,-1)
                if PrF < .0001:
                    output.append((Fenter,k))
                continue
        break

    output.sort(reverse=True)
    return array([x[1] for x in output])
