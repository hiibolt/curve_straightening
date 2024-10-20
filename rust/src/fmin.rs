fn fmin (
    ax: f64,
    bx: f64,
    tol: f64
) -> f64 {
    // C is the squared inverse of the golden ratio.
    let c = (3.0 - (5.0f64).sqrt()) / 2.0;

    // Initialize.
    let a = ax;
    let b = bx;
    let v = a + c * (b - a);
    let w = v;
    let x = v;
    let e = 0.0;

    let mut fmin = x;
    let mut iflg = 0;

    // EPS is approximately the square root of the smallest
    // positive floating-point number such that 1.0 + EPS != 1.0.
    let eps = std::f64::EPSILON.sqrt();

    // Reason for todo: phi is not implemented.
    todo!();

    return 0.0;
}