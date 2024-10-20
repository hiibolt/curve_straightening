use crate::feval::{ feval, FevalResult };
use crate::initc::initc;

/* 
C***********************************************************
C
C                                                  From GEO2
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/19/03
C
C   Given a polygonal curve defined by a sequence of verti-
C ces c(i), i = 0 to m, along with a weight w associated
C with the constraint penalty, this function returns the
C following:
C
C       ds(i) = |c(i)-c(i-1)|,
C
C       D1c(i) = [c(i)-c(i-1)]/ds(i),
C
C       D2c(i) = [D1c(i+1)-D1c(i)]/([ds(i)+ds(i+1)]/2)
C
C for i = 1 to m, where ds(m+1) = ds(1) and D1c(m+1) =
C D1c(1), and the total geodesic curvature
C
C       E(c) = 0.5*Sum[ (|u(i)|**2)*da(i) ]
C
C where
C
C       da(i) = [ds(i)+ds(i+1)],
C
C       u(i) = [I - fn(i)*fn(i)**T]*D2c(i)
C            = D2c(i) - <fn(i),D2c(i)>*fn(i),
C
C       fn(i) = grad(f)(c(i))/|grad(f)(c(i))|,
C
C and Sum[_] denotes the sum over i = 1 to m.  This function
C also returns
C
C       FNS = 0.5*Sum[ F(i)**2 ]
C
C and
C
C       phi(c) = E(c) + w*FNS
C
C where F(i) = f(c(i)).
C
C Refer to the comments at the beginning of the main program
C for more details.
C
C   Subroutine FEVAL is called to compute F(i) = f(c(i)) and
C the gradients of f at c(i) for i = 1 to m.
C
C On input:
C
C       M = Number of vertices and curve segments.  M >= 3.
C
C       C = Array dimensioned 3 by 0:M containing the
C           vertices with the Cartesian coordinates of
C           c(j) in column j.  Adjacent vertices c(i-1)
C           and c(i) must be distinct for i = 1 to m.
C           It is assumed that c(0) = c(m).
C
C       W = Weight defining the penalty term.  W >= 0.
C
C The above parameters are not altered by this function.
C
C       DS = Array of length at least M.
C
C       D1C,D2C,U = Arrays dimensioned 3 by M.
C
C       F = Array of length at least M.
C
C       G,HU = Arrays dimensioned 3 by M.
C
C On output:
C
C       DS,D1C,D2C,U = Arrays containing the terms defined
C                      above:  arc lengths, unit tangent
C                      vectors, curvature vectors, and cur-
C                      vature vectors projected onto tangent
C                      planes.
C
C       F,G = Arrays containing function values f(c(i)) and
C             gradient vectors grad(f)(c(i)), respectively,
C             computed by Subroutine FEVAL.
C
C       HU = Array containing terms needed for the gradient
C            of phi:  s*H(i)*u(i), for i = 1 to m, where
C            s = -da(i)*<fn(i),D2c(i)>/|G( ,i)| and H(i) is
C            the Hessian of f at c(i).
C
C       CL = Total curve length:  Sum[ ds(i) ].
C
C       E = Value of the functional representing total
C           geodesic curvature.
C
C       FNS = Half the sum of squared components of F.
C
C       PHI = E(c) + w*FNS.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if M or W is outside its valid range
C                     on input.  Output parameters are not
C                     altered in this case.
C             IER = 2 if c(i-1) = c(i) for some i (1 to m).
C             IER = 3 if grad(f)(c(i)) = 0 for some i
C                     (1 to m).
C
C Module required by PHI:  FEVAL
C
C Intrinsic function called by PHI:  SQRT
C
C***********************************************************
*/
#[derive(Debug)]
pub enum PhiError {
    InvalidMW,
    DuplicateVertices,
    ZeroGradient,
}
pub struct PhiResult {
    pub cl: f64,
    pub e: f64,
    pub fns: f64,
    pub phi: f64,
}
fn phi (
    m: usize,
    c: &Vec<Vec<f64>>,
    w: f64,

    ds: &mut Vec<f64>,
    d1c: &mut Vec<Vec<f64>>,
    d2c: &mut Vec<Vec<f64>>,
    u: &mut Vec<Vec<f64>>,
    f: &mut Vec<f64>,
    g: &mut Vec<Vec<f64>>,
    hu: &mut Vec<Vec<f64>>,
) -> Result<PhiResult, PhiError> {
    // Check for invalid input
    if m < 3 || w < 0.0 {
        return Err(PhiError::InvalidMW);
    }

    // Compute ds(i) and D1c(i) for i = 1 to m, where
    //  Sum[ds(i)] is also computed and stored in CL.
    let mut cl = 0.0f64;
    for i in 1..=m {
        // Compute arc lengths.
        ds[i - 1] = (c[0][i] - c[0][i - 1]).powi(2) + 
                    (c[1][i] - c[1][i - 1]).powi(2) + 
                    (c[2][i] - c[2][i - 1]).powi(2);
        
        // Check for duplicate vertices.
        if ds[i - 1] == 0.0 {
            return Err(PhiError::DuplicateVertices);
        }

        ds[i - 1] = ds[i - 1].sqrt();
        cl += ds[i - 1];

        d1c[0][i - 1] = (c[0][i] - c[0][i - 1]) / ds[i - 1];
        d1c[1][i - 1] = (c[1][i] - c[1][i - 1]) / ds[i - 1];
        d1c[2][i - 1] = (c[2][i] - c[2][i - 1]) / ds[i - 1];
    }

    /*
    C Compute D2c(i),
    C         F(i) = f(c(i)),
    C         G( ,i) = grad(f)(c(i)),
    C         u(i) = D2c(i) - <fn(i),D2c(i)>*fn(i), and
    C         HU( ,i) = s*H(i)*u(i),
    C
    C for i = 1 to m-1, where s = -da(i)*<fn(i),D2c(i)>/|G( ,i)|
    C and H(i) is the Hessian of f at c(i).
    C
    C Sum[ (|u(i)|**2)*da(i) ] is accumulated in CSUM, and
    C Sum[ F(i)**2 ] is accumulated in FSUM.
    */
    let mut csum = 0.0f64;
    let mut fsum = 0.0f64;
    let mut count = 0;
    for i in 1..=(m-1) {
        let da = ds[i - 1] + ds[i];
        let mut s = 2.0f64 / da;

        d2c[0][i - 1] = s * (d1c[0][i] - d1c[0][i - 1]);
        d2c[1][i - 1] = s * (d1c[1][i] - d1c[1][i - 1]);
        d2c[2][i - 1] = s * (d1c[2][i] - d1c[2][i - 1]);

        let feval_result = feval(c[0][i], c[1][i], c[2][i], 2);
        match feval_result {
            FevalResult::SecondOrder { 
                f: fval, fx, fy, fz, fxx, fxy, fxz, fyy, fyz, fzz
            } => {
                f[i - 1] = fval;
                g[0][i - 1] = fx;
                g[1][i - 1] = fy;
                g[2][i - 1] = fz;

                fsum += fval.powi(2);

                s = g[0][i - 1].powi(2) + g[1][i - 1].powi(2) + g[2][i - 1].powi(2);

                if s == 0.0 {
                    return Err(PhiError::ZeroGradient);
                }

                s = (g[0][i - 1] * d2c[0][i - 1] + 
                     g[1][i - 1] * d2c[1][i - 1] + 
                     g[2][i - 1] * d2c[2][i - 1]) / s;
                
                u[0][i - 1] = d2c[0][i - 1] - s * g[0][i - 1];
                u[1][i - 1] = d2c[1][i - 1] - s * g[1][i - 1];
                u[2][i - 1] = d2c[2][i - 1] - s * g[2][i - 1];

                csum += da * (u[0][i - 1].powi(2) + 
                              u[1][i - 1].powi(2) + 
                              u[2][i - 1].powi(2));
                count += 1;

                s = -da * s;

                hu[0][i - 1] = s * fxx * u[0][i - 1] + 
                           s * fxy * u[1][i - 1] + 
                           s * fxz * u[2][i - 1];
                hu[1][i - 1] = s * fxy * u[0][i - 1] +
                           s * fyy * u[1][i - 1] +
                           s * fyz * u[2][i - 1];
                hu[2][i - 1] = s * fxz * u[0][i - 1] +
                           s * fyz * u[1][i - 1] +
                           s * fzz * u[2][i - 1];
            },
            _ => panic!("Unexpected result type")
        }
    }

    // Compute D2c(m), F(m), G( ,m), u(m), and HU( ,m), and add
    // the final contribution to FSUM and CSUM.
    let da = ds[m - 1] + ds[0];
    let mut s = 2.0f64 / da;

    d2c[0][m - 1] = s * (d1c[0][0] - d1c[0][m - 1]);
    d2c[1][m - 1] = s * (d1c[1][0] - d1c[1][m - 1]);
    d2c[2][m - 1] = s * (d1c[2][0] - d1c[2][m - 1]);

    let feval_result = feval(c[0][m], c[1][m], c[2][m], 2);
    match feval_result {
        FevalResult::SecondOrder { 
            f: fval, fx, fy, fz, fxx, fxy, fxz, fyy, fyz, fzz
        } => {
            f[m - 1] = fval;
            g[0][m - 1] = fx;
            g[1][m - 1] = fy;
            g[2][m - 1] = fz;

            fsum += fval.powi(2);

            s = g[0][m - 1].powi(2) + 
                g[1][m - 1].powi(2) + 
                g[2][m - 1].powi(2);
            
            if s == 0.0 {
                return Err(PhiError::ZeroGradient);
            }

            s = (g[0][m - 1] * d2c[0][m - 1] + 
                 g[1][m - 1] * d2c[1][m - 1] + 
                 g[2][m - 1] * d2c[2][m - 1]) / s;

            u[0][m - 1] = d2c[0][m - 1 ] - s * g[0][m - 1];
            u[1][m - 1] = d2c[1][m - 1 ] - s * g[1][m - 1];
            u[2][m - 1] = d2c[2][m - 1 ] - s * g[2][m - 1];

            csum += da * (u[0][m - 1].powi(2) + 
                          u[1][m - 1].powi(2) + 
                          u[2][m - 1].powi(2));
            s = -da * s;

            hu[0][m - 1] = s * fxx * u[0][m - 1] + 
                           s * fxy * u[1][m - 1] + 
                           s * fxz * u[2][m - 1];
            hu[1][m - 1] = s * fxy * u[0][m - 1] +
                           s * fyy * u[1][m - 1] +
                           s * fyz * u[2][m - 1];
            hu[2][m - 1] = s * fxz * u[0][m - 1] +
                           s * fyz * u[1][m - 1] +
                           s * fzz * u[2][m - 1];
        },
        _ => panic!("Unexpected result type")
    }

    // Compute E, FNS, and PHI.
    let e = 0.5f64 * csum;
    let fns = 0.5 * fsum;
    let phi = e + w * fns;

    Ok(PhiResult {
        cl,
        e,
        fns,
        phi
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let m = 400;
        let mut c: Vec<Vec<f64>> = vec![
            vec![0.0; m + 1],
            vec![0.0; m + 1],
            vec![0.0; m + 1]
        ];

        initc(m, &mut c);

        let result = phi(
            m, &c, 0.0, 
            &mut vec![0.0; m], 
            &mut vec![vec![0.0; m]; 3], 
            &mut vec![vec![0.0; m]; 3], 
            &mut vec![vec![0.0; m]; 3], 
            &mut vec![0.0; m], 
            &mut vec![vec![0.0; m]; 3], 
            &mut vec![vec![0.0; m]; 3]
        ).unwrap();

        assert_eq!(result.cl,  5.1506186431346039);
        assert_eq!(result.e,   22.839762081427391);
        assert_eq!(result.fns, 99.109515912078862);
        assert_eq!(result.phi, 22.839762081427391);
    }
}