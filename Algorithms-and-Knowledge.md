# Algorithms and Knowledge

## Normalization

Giving a group of cells:
$$ \left\{ i \mid i\in \mathbb{C}_g \right\} $$

Reference cells:
$$ \left\{ r \mid r\in \mathbb{C}_{gr} \,\text{and}\, \mathbb{C}_{gr}\in\mathbb{C}_g \right\} $$

Log expression difference:
$$ r_{ir\_j} = \log({e_{ij}}) - \log({e_{rj}}) $$

Scaling factor:
$$ s_{ir} = \sum_{j,j\in\mathbb{S}_J}= \left. r_{ir_j} / \right.  \left|\mathbb{S}_J\right| $$ where $\left|\mathbb{S}_J\right|$ is the number of genes

Relate different reference cells:
$$ \left\{ \hat{a}_r \mid r\in\mathbb{C}_{gr} \right\} = \mathop{\arg\max}\limits_{a_r,r\in\mathbb{C}_{gr}} \sum_{i\in\mathbb{C}_{gr}} \mathop{Var}\limits_{r\in\mathbb{C}_{gr}}[\log(s_{ir}) - \log(a_r)] $$

Estimated scaling factors:
$$ s_i = \mathop{median}\limits_{r\in\mathbb{C}_{gr}} \left( \frac{s_{ir}}{\hat{a}_r} \right) $$

Expression of gene *j* in group *g* :
$$ e_{gi} = \sum_{i\in\mathbb{C}_g} e'_{ij}$$
where $e'_{ij}$ is the normalized expression level.

## Kernel density estimation

The kernel density estimate is defined to be
$$ \hat{f_{\mathrm{H}}}(\mathrm{x}) = \frac{1}{n}\sum\limits^n\limits_{i=1}K_{\mathrm{H}}(\mathrm{x} - \mathrm{x}_i) $$
where
* $\mathrm{x}=(x_1,x_2,...,x_d)^T,x_i=(x_{i1},x_{i2},..x_{id})^T,i=1,2,...,n$ are *d*-vectors;
* $H$ is the bandwith (or smooting) *d×d* matrix which is symmetric and positive definite;
* *K* is the kernel function which is a symmetric multivariate density;
* $K_{\mathrm{H}}(\mathrm{x})=|\mathrm{H}|^{-1/2}K(\mathrm{H}^{-1/2}\mathrm{x})$

