# SPDX-License-Identifier: Apache-2.0
"""Validation metrics and site-level pooling, stdlib only.

- observed/expected ratio, Brier score, C-statistic (binary outcome, ties 1/2)
- logistic calibration slope and intercept (Newton-Raphson) and
  calibration-in-the-large (intercept with logit(p) as offset)
- time-to-event observed risk by a horizon: Aalen-Johansen cumulative incidence
  (other causes compete) and 1 - Kaplan-Meier (other causes censored), and the
  O/E ratio built on either (new in 0.2.1)
- DerSimonian-Laird random-effects pooling with I^2 and the Higgins-Thompson-
  Spiegelhalter (2009) 95% prediction interval for a NEW site (t, k-2 df)
- random-effects meta-regression on one site-level moderator: method-of-moments
  tau^2, then Knapp-Hartung standard errors with a t(k-2) test

Research code only. Information only, not medical advice, not a medical device.
"""
import math

# Two-sided 95% t quantiles t_{0.975, df}, df = 1..30
_T975 = [12.706205, 4.302653, 3.182446, 2.776445, 2.570582, 2.446912, 2.364624,
         2.306004, 2.262157, 2.228139, 2.200985, 2.178813, 2.160369, 2.144787,
         2.131450, 2.119905, 2.109816, 2.100922, 2.093024, 2.085963, 2.079614,
         2.073873, 2.068658, 2.063899, 2.059539, 2.055529, 2.051831, 2.048407,
         2.045230, 2.042272]
_Z975 = 1.959964


def t975(df):
    """t_{0.975, df}: table for df <= 30, Cornish-Fisher expansion above."""
    if df < 1:
        raise ValueError("df must be >= 1")
    if df <= 30:
        return _T975[int(df) - 1]
    z = _Z975
    return (z + (z ** 3 + z) / (4.0 * df)
            + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96.0 * df * df))


def logit(p):
    return math.log(p / (1.0 - p))


def expit(x):
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def observed_expected(y, p):
    e = sum(p)
    if e <= 0:
        raise ValueError("sum of predicted risks must be > 0")
    return sum(y) / e


def brier(y, p):
    if not y:
        raise ValueError("empty")
    return sum((pi - yi) ** 2 for yi, pi in zip(y, p)) / len(y)


def c_statistic(y, p):
    """Probability that a random event has a higher predicted risk than a random
    non-event (ties count 1/2). O(n log n) via ranks."""
    pairs = sorted(zip(p, y))
    n = len(pairs)
    n1 = sum(1 for _, yi in pairs if yi)
    n0 = n - n1
    if n1 == 0 or n0 == 0:
        return None
    rank_sum = 0.0
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        rank_sum += avg_rank * sum(1 for k in range(i, j + 1) if pairs[k][1])
        i = j + 1
    return (rank_sum - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def calibration_slope(y, p, max_iter=100, tol=1e-12):
    """Logistic regression of y on logit(p). Returns
    (intercept, slope, se_intercept, se_slope)."""
    x = [logit(pi) for pi in p]
    a, b = 0.0, 1.0
    for _ in range(max_iter):
        g0 = g1 = h00 = h01 = h11 = 0.0
        for xi, yi in zip(x, y):
            m = expit(a + b * xi)
            w = m * (1.0 - m)
            g0 += yi - m
            g1 += (yi - m) * xi
            h00 += w
            h01 += w * xi
            h11 += w * xi * xi
        det = h00 * h11 - h01 * h01
        if det <= 0:
            raise ValueError("calibration model not identifiable")
        da = (h11 * g0 - h01 * g1) / det
        db = (h00 * g1 - h01 * g0) / det
        a += da
        b += db
        if abs(da) < tol and abs(db) < tol:
            break
    h00 = h01 = h11 = 0.0
    for xi in x:
        m = expit(a + b * xi)
        w = m * (1.0 - m)
        h00 += w
        h01 += w * xi
        h11 += w * xi * xi
    det = h00 * h11 - h01 * h01
    return a, b, math.sqrt(h11 / det), math.sqrt(h00 / det)


def calibration_in_the_large(y, p, max_iter=100, tol=1e-12):
    """Intercept a in logit P(y=1) = a + logit(p) (slope fixed at 1).
    Returns (a, se_a). a > 0 means the model under-predicts."""
    x = [logit(pi) for pi in p]
    a = 0.0
    for _ in range(max_iter):
        g = h = 0.0
        for xi, yi in zip(x, y):
            m = expit(a + xi)
            g += yi - m
            h += m * (1.0 - m)
        step = g / h
        a += step
        if abs(step) < tol:
            break
    h = sum(expit(a + xi) * (1.0 - expit(a + xi)) for xi in x)
    return a, math.sqrt(1.0 / h)


# ------------------------------------------------------------ time to event
def cumulative_incidence(times, events, horizon, cause=1):
    """Aalen-Johansen cumulative incidence of `cause` by time `horizon`.

    events: 0 = censored, any other value = an event of that type; every
    non-zero type other than `cause` is a competing event. Someone censored at
    time t is still at risk at t (events before censorings at tied times).
    Returns (cumulative incidence, number of `cause` events at or before the
    horizon)."""
    data = sorted(zip(times, events))
    if not data:
        raise ValueError("empty")
    at_risk = len(data)
    surv = 1.0          # all-cause event-free probability just before t
    cif = 0.0
    n_cause = 0
    i = 0
    while i < len(data):
        t = data[i][0]
        if t > horizon:
            break
        j = i
        d_cause = d_any = 0
        while j < len(data) and data[j][0] == t:
            e = data[j][1]
            if e == cause:
                d_cause += 1
            if e != 0:
                d_any += 1
            j += 1
        if d_any:
            cif += surv * d_cause / at_risk
            surv *= 1.0 - d_any / at_risk
        n_cause += d_cause
        at_risk -= j - i
        i = j
    return cif, n_cause


def km_risk(times, events, horizon, cause=1):
    """1 - Kaplan-Meier survival from `cause` by `horizon`, every other event
    type treated as censoring. Larger than the Aalen-Johansen cumulative
    incidence whenever a competing event happens before a `cause` event.
    Returns (risk, number of `cause` events at or before the horizon)."""
    return cumulative_incidence(times, [1 if e == cause else 0 for e in events],
                                horizon, cause=1)


def oe_time_to_event(times, events, p, horizon, cause=1, competing=True):
    """O/E by a horizon for everyone with a prediction, including people whose
    follow-up ended early. O = Aalen-Johansen cumulative incidence (competing=True)
    or 1 - Kaplan-Meier (competing=False); E = mean predicted risk.
    se_log_oe = 1/sqrt(events), a Poisson approximation (None with 0 events).
    Returns a dict."""
    if not p or len(p) != len(times) or len(times) != len(events):
        raise ValueError("need matching, non-empty times, events and predictions")
    if competing:
        risk, n_ev = cumulative_incidence(times, events, horizon, cause)
    else:
        risk, n_ev = km_risk(times, events, horizon, cause)
    expected = sum(p) / len(p)
    if expected <= 0:
        raise ValueError("mean predicted risk must be > 0")
    return {"n": len(p), "events": n_ev, "observed_risk": risk, "expected_risk": expected,
            "oe": risk / expected,
            "se_log_oe": (1.0 / math.sqrt(n_ev)) if n_ev > 0 else None}


# ------------------------------------------------------------ pooling
def dersimonian_laird(estimates, variances):
    """Random-effects pooling. Returns a dict with pooled, se, tau2, q, i2,
    ci (95%, normal) and pi (95% prediction interval for a new site, t with
    k-2 df; None when k < 3)."""
    k = len(estimates)
    if k < 2 or k != len(variances):
        raise ValueError("need >= 2 estimates with matching variances")
    w = [1.0 / v for v in variances]
    sw = sum(w)
    fixed = sum(wi * yi for wi, yi in zip(w, estimates)) / sw
    q = sum(wi * (yi - fixed) ** 2 for wi, yi in zip(w, estimates))
    c = sw - sum(wi * wi for wi in w) / sw
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    ws = [1.0 / (v + tau2) for v in variances]
    pooled = sum(wi * yi for wi, yi in zip(ws, estimates)) / sum(ws)
    se = math.sqrt(1.0 / sum(ws))
    ci = (pooled - _Z975 * se, pooled + _Z975 * se)
    pi = None
    if k >= 3:
        half = t975(k - 2) * math.sqrt(tau2 + se * se)
        pi = (pooled - half, pooled + half)
    return {"k": k, "pooled": pooled, "se": se, "tau2": tau2, "q": q, "i2": i2,
            "ci": ci, "pi": pi}


def _wls(x, y, w):
    """Weighted least squares y = b0 + b1 x. Returns b0, b1 and the inverse
    of X'WX as (m00, m01, m11)."""
    s0 = sum(w)
    s1 = sum(wi * xi for wi, xi in zip(w, x))
    s2 = sum(wi * xi * xi for wi, xi in zip(w, x))
    t0 = sum(wi * yi for wi, yi in zip(w, y))
    t1 = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
    det = s0 * s2 - s1 * s1
    if det <= 0:
        raise ValueError("moderator has no spread")
    b0 = (s2 * t0 - s1 * t1) / det
    b1 = (s0 * t1 - s1 * t0) / det
    return b0, b1, (s2 / det, -s1 / det, s0 / det)


def meta_regression(estimates, variances, moderator):
    """Random-effects meta-regression of site estimates on one moderator.
    tau^2: method of moments (DerSimonian-Laird type) from the fixed-effect
    residual Q. Inference: Knapp-Hartung, t with k-2 df.
    Returns a dict with b0, b1, se_b1, t, df, ci_b1, reject (two-sided 5%),
    tau2, q_res."""
    k = len(estimates)
    if k < 3:
        raise ValueError("need >= 3 sites")
    w = [1.0 / v for v in variances]
    b0, b1, (m00, m01, m11) = _wls(moderator, estimates, w)
    q_res = sum(wi * (yi - b0 - b1 * xi) ** 2
                for wi, xi, yi in zip(w, moderator, estimates))
    # tr[(X'WX)^-1 X'W^2X]
    u0 = sum(wi * wi for wi in w)
    u1 = sum(wi * wi * xi for wi, xi in zip(w, moderator))
    u2 = sum(wi * wi * xi * xi for wi, xi in zip(w, moderator))
    tr = m00 * u0 + 2.0 * m01 * u1 + m11 * u2
    denom = sum(w) - tr
    tau2 = max(0.0, (q_res - (k - 2)) / denom) if denom > 0 else 0.0
    ws = [1.0 / (v + tau2) for v in variances]
    b0, b1, (m00, m01, m11) = _wls(moderator, estimates, ws)
    s2 = sum(wi * (yi - b0 - b1 * xi) ** 2
             for wi, xi, yi in zip(ws, moderator, estimates)) / (k - 2)
    se_b1 = math.sqrt(s2 * m11)
    df = k - 2
    tq = t975(df)
    t = b1 / se_b1 if se_b1 > 0 else (0.0 if b1 == 0 else math.copysign(math.inf, b1))
    return {"k": k, "b0": b0, "b1": b1, "se_b1": se_b1, "t": t, "df": df,
            "ci_b1": (b1 - tq * se_b1, b1 + tq * se_b1),
            "reject": abs(t) > tq, "tau2": tau2, "q_res": q_res}
