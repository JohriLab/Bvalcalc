import scipy.stats as st

def gammaDFE_to_discretized(mean: float, shape: float):
    """
    Discretise a Gamma distribution into four probability bins:
        f0: 0‑1
        f1: 1‑10
        f2: 10‑100
        f3: >100

    Parameters
    ----------
    mean  : float
        Mean (μ) of the Gamma distribution.
    shape : float
        Shape parameter (k or α).

    Returns
    -------
    tuple[float, float, float, float]
        (f0, f1, f2, f3) – the probability mass in each interval.
    """
    if mean <= 0 or shape <= 0:
        raise ValueError("Both mean and shape must be positive.")

    theta = mean / shape               # scale parameter
    dist  = st.gamma(a=shape, scale=theta)

    # cumulative probabilities at the cut‑points
    c1   = dist.cdf(1.0)
    c10  = dist.cdf(10.0)
    c100 = dist.cdf(100.0)

    f0 = c1                      # (0,1]
    f1 = c10  - c1               # (1,10]
    f2 = c100 - c10              # (10,100]
    f3 = 1.0   - c100            # >100

    return f0, f1, f2, f3
