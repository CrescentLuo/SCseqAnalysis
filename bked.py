"""
An implementation of the kde bandwidth selection method outlined in:

Reference: Z. I. Botev, J. F. Grotowski and D. P. Kroese
            "KERNEL DENSITY ESTIMATION VIA DIFFUSION" ,Submitted to the
            Annals of Statistics, 2009

Based on the implementation in Matlab by Zdravko Botev.
"""

import scipy
import numpy as np
from numpy.matlib import repmat
from scipy.fftpack import dct, idct


def kde2d(data, grid_n):
    grid_n = 2** 8 if N is None else int(2**np.ceil(sci.log2(N)))
    data_n = data.shape[0]
    data_max = np.max(data, axis=0)
    data_min = np.min(data, axis=0)
    data_range = data_max - data_min
    data_max_xy = data_max + data_range/4
    data_min_xy = data_min - data_range/4
    scaling = data_max_xy - data_min_xy

    transformed_data=np.divide(data - repmat(data_min_xy, data_n, 1), repmat(scaling, data_n, 1))
    #bin the data uniformly using regular grid;
    init_data = np.histogram2d(transformed_data[0], transformed_data[1], grid_n)
    # discrete cosine transform of initial data
    SqrtDCTData = dct(init_data) ** 2
    # now compute the optimal bandwidth^2
    I = (np.arange(0, grid_n) **2).reshape(1,-1)
    SqrtDCTData = a**2
    t_star = scipy.optimize.brentq(fixed_point, 0, 0.1, 
                                        args=(M, I, SqDCTData))
def fixed_point(t, M, I, a2):
        l=7
        I = sci.float128(I)
        M = sci.float128(M)
        a2 = sci.float128(a2)
        f = 2*sci.pi**(2*l)*sci.sum(I**l*a2*sci.exp(-I*sci.pi**2*t))
        for s in range(l, 1, -1):
            K0 = sci.prod(xrange(1, 2*s, 2))/sci.sqrt(2*sci.pi)
            const = (1 + (1/2)**(s + 1/2))/3
            time=(2*const*K0/M/f)**(2/(3+2*s))
            f=2*sci.pi**(2*s)*sci.sum(I**s*a2*sci.exp(-I*sci.pi**2*time))
        return t-(2*M*sci.sqrt(sci.pi)*f)**(-2/5)
    t_star = 
    t_star=fzero( @(t)(t-evolve(t)),[0,0.1]);

    p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
    t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    # smooth the discrete cosine transform of initial data using t_star
    a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
    # now apply the inverse discrete cosine transform
    if nargout>1
        density=idct2d(a_t)*(numel(a_t)/prod(scaling));
        [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
    end
    bandwidth=sqrt([t_x,t_y]).*scaling; 
    end

    def evolve(t):
        Sum_func = func([0, 2], t) + func([2, 0], t) + 2* func([1, 1], t)
        time = (2 * np.pi * grid_n * Sum_func)^(-1/3)
    ########################################
    function  [out,time]=evolve(t)
    global N
    Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
    time=(2*pi*N*Sum_func)^(-1/3);
    out=(t-time)/time;
    end
    ########################################
    function out=func(s,t)
    global N
    if sum(s)<=4
        Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
        time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
        out=psi(s,time);
    else
        out=psi(s,t);
    end

    end
    ########################################
    def psi(s, time)：

    global I, SqrtDCTData
    # s is a vector
    w = np.exp(-I * np.pi ** 2 * time) * (0.5 * np.ones((1, len(I)-1)))
    w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
    wx=w.*(I.^s(1));
    wy=w.*(I.^s(2));
    out=(-1)^sum(s)*(wy*SqrtDCTData*wx')*pi^(2*sum(s));
    return out
    ########################################
    function out=K(s)
    out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
    end
    ########################################
    function binned_data=ndhist(data,M)
    # this function computes the histogram
    # of an n-dimensional data set;
    # 'data' is nrows by n columns
    # M is the number of bins used in each dimension
    # so that 'binned_data' is a hypercube with
    # size length equal to M;
    [nrows,ncols]=size(data);
    bins=zeros(nrows,ncols);
    for i=1:ncols
        [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
        bins(:,i) = min(bins(:,i),M);
    end
    # Combine the  vectors of 1D bin counts into a grid of nD bin
    # counts.
    binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
    end
