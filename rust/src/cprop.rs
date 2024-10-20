#[derive(Debug)]
pub struct CPROPResult {
    pub cnrm2: f64,
    pub cmax: f64,
}

/**
C***********************************************************
C
C                                                  From GEO1
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/19/03
C
C   Given a polygonal approximation to a closed curve c in
C a level set {(x,y,z):  f(x,y,z) = 0}, along with segment
C lengths ds(i) and discrete approximations u(i) to vertex
C values of geodesic curvature Kg, this subroutine computes
C vertex values |u(i)| and norms of the total geodesic
C curvature.
C
C On input:
C
C       m (num_segments) = Number of vertices and polygonal curve segments
C           in the approximation of c.  m >= 1.
C
C       ds (segment_lengths) = Array of length m containing segment lengths:
C            ds(i) = |c(i)-c(i-1)| for i = 1 to m.
C
C       u = Array dimensioned 3 by m containing vertex
C           values of [I - fn(i)*fn(i)**T]*D2c(i), where
C           fn(i) is the unit normal to the surface and
C           D2c(i) is an approximation to c'' at c(i).
C
C The above parameters are not altered by this routine.
C
C       s,crv = Arrays of length at least m.
C
C On output:
C
C       s = Cumulative segment lengths:  s(j) = Sum[ds(i)],
C           where the sum is over i = 1 to j for j = 1 to m.
C           For the natural parameterization by arc length,
C           s(j) is the parameter value associated with ver-
C           tex j, s(j)-ds(j)/2 is the parameter value
C           associated with midpoint j, and s(m) is the
C           total curve length.
C
C       crv = Vertex values of the geodesic curvature
C             |Kg(s)| = |u(i)|.
C
C       cnrm2 = Approximation to the integral with respect
C               to arc length s of Kg(s)**2.
C
C       cmax = Approximation to the max-norm of |Kg(s)|.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if m < 1.
*/
pub fn cprop (
    m: usize,
    ds: &Vec<f64>,
    u: &Vec<Vec<f64>>,

    s: &mut Vec<f64>,
    crv: &mut Vec<f64>,
) -> Result<CPROPResult, ()> {
    if m < 1 {
        return Err(());
    }

    // Compute cumulative arc length values s(j), j = 1 to m.
    s[0] = ds[0];
    for i in 1..m {
        s[i] = s[i-1] + ds[i];
    }

    // Compute vertex geodesic curvature values Kg(s).
    for i in 0..m {
        crv[i] = (u[i][0].powi(2) + u[i][1].powi(2) + u[i][2].powi(2)).sqrt();
    }

    // Accumulate norms using the composite midpoint rule for
    //  cnrm2.  cj is the midpoint value of Kg(s).
    let mut cnrm2: f64 = 0.0;
    let mut cmax: f64 = 0.0;
    let mut cj: f64 = (crv[m-1] + crv[0]) / 2.0;
    for i in 0..m {
        cnrm2 += ds[i] * cj.powi(2);
        cmax = cmax.max(cj);

        if i < m-1 {
            cj = (crv[i] + crv[i+1]) / 2.0f64;
        }
    }

    Ok(CPROPResult {
        cnrm2,
        cmax,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let m: usize = 4;
        let ds: Vec<f64> = vec![1.0, 1.5, 2.0, 2.5];
        let u: Vec<Vec<f64>> = vec![
            vec![0.1, 0.2, 0.3],
            vec![0.4, 0.5, 0.6],
            vec![0.7, 0.8, 0.9],
            vec![1.0, 1.1, 1.2]
        ];

        let mut s: Vec<f64> = vec![0.0; m + 1];
        let mut crv: Vec<f64> = vec![0.0; m + 1];

        let result = cprop(m, &ds, &u, &mut s, &mut crv).unwrap();

        assert_eq!(crv[0], 0.37416573867739417);
        assert_eq!(crv[1], 0.87749643873921224);
        assert_eq!(crv[2], 1.3928388277184121);
        assert_eq!(crv[3], 1.9104973174542801);

        assert_eq!(result.cnrm2, 11.289647813467305);
        assert_eq!(result.cmax, 1.6516680725863462);
    }
}