/// Enum to specify the different result levels
pub enum FevalResult {
    ValueOnly(f64),
    FirstOrder {
        f: f64,
        fx: f64,
        fy: f64,
        fz: f64,
    },
    SecondOrder {
        f: f64,
        fx: f64,
        fy: f64,
        fz: f64,
        fxx: f64,
        fxy: f64,
        fxz: f64,
        fyy: f64,
        fyz: f64,
        fzz: f64,
    },
}

/// Evaluates the function and its partial derivatives
///
/// # Arguments
///
/// * `x`, `y`, `z` - The evaluation points.
/// * `level` - The level of derivatives to return (ValueOnly, FirstOrder, or SecondOrder).
///
/// # Returns
///
/// Returns a `FevalResult` enum with computed values based on `level`.
pub fn feval(x: f64, y: f64, z: f64, level: u8) -> FevalResult {
    // Compute the function value
    let f = 2.0f64*x.powi(2) + 3.0f64*y.powi(2) + 5.0f64*z.powi(2) + x.powi(3)*z + y*z.powi(3) - 1.0f64;

    // Compute the first partial derivatives
    let (fx, fy, fz) = if level > 0 {
        (
            4.0f64 * x + 3.0f64 * x.powi(2) * z,
            6.0f64 * y + z.powi(3),
            10.0f64 * z + x.powi(3) + 3.0f64 * y * z.powi(2),
        )
    } else {
        (0.0, 0.0, 0.0)
    };

    // Compute the second partial derivatives
    let (fxx, fxy, fxz, fyy, fyz, fzz) = if level > 1 {
        (
            4.0f64 + 6.0f64 * x * z,
            0.0,
            3.0f64 * x.powi(2),
            6.0,
            3.0f64 * z.powi(2),
            10.0f64 + 6.0f64 * y * z,
        )
    } else {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    };

    // Return the appropriate result based on the level
    if level == 0 {
        FevalResult::ValueOnly(f)
    } else if level == 1 {
        FevalResult::FirstOrder { f, fx, fy, fz }
    } else {
        FevalResult::SecondOrder {
            f,
            fx,
            fy,
            fz,
            fxx,
            fxy,
            fxz,
            fyy,
            fyz,
            fzz,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_level_one() {
        let x: f64 = 1.0;
        let y: f64 = 1.0;
        let z: f64 = 1.0;
        let level: u8 = 0;

        let result = feval(x, y, z, level);

        match result {
            FevalResult::ValueOnly(f) => {
                assert_eq!(f, 11.000000000000000);
            }
            _ => panic!("Unexpected result type"),
        }
    }

    #[test]
    fn basic_level_two() {
        let x: f64 = 1.0;
        let y: f64 = 1.0;
        let z: f64 = 1.0;
        let level: u8 = 1;

        let result = feval(x, y, z, level);

        match result {
            FevalResult::FirstOrder { f, fx, fy, fz } => {
                assert_eq!(f,  11.000000000000000);
                assert_eq!(fx, 7.0000000000000000);
                assert_eq!(fy, 7.0000000000000000);
                assert_eq!(fz, 14.000000000000000);
            }
            _ => panic!("Unexpected result type"),
        }
    }

    #[test]
    fn basic_level_three() {
        let x: f64 = 1.0;
        let y: f64 = 1.0;
        let z: f64 = 1.0;
        let level: u8 = 2;

        let result = feval(x, y, z, level);

        match result {
            FevalResult::SecondOrder { f, fx, fy, fz, fxx, fxy, fxz, fyy, fyz, fzz } => {
                assert_eq!(f,   11.000000000000000);
                assert_eq!(fx,  7.0000000000000000);
                assert_eq!(fy,  7.0000000000000000);
                assert_eq!(fz,  14.000000000000000);
                assert_eq!(fxx, 10.000000000000000);
                assert_eq!(fxy, 0.0000000000000000);
                assert_eq!(fxz, 3.0000000000000000);
                assert_eq!(fyy, 6.0000000000000000);
                assert_eq!(fyz, 3.0000000000000000);
                assert_eq!(fzz, 16.000000000000000);
            }
            _ => panic!("Unexpected result type"),
        }
    }
}