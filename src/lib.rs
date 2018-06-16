extern crate rand;

// randvariate - module that provides functions to generate random variates
mod randvariate {
    use rand::distributions::{Uniform, Distribution};
    use rand::{ThreadRng, thread_rng};

    const CON:f64 = 57.56462733;
    const DELTAL:f64 = 0.0078;
    const DELTAU:f64 = 0.0034;
    const SCALE:f64 = 1.0e25;


    /// h2pec - Hypergeometric Random Variate Generator
    /// given parameters: kk (number draws) nn1 (number white balls)
    /// nn2 (number black balls) and rng (random number generator)
    /// hgd will return the number of white balls drawn
    ///
    /// This function is based on the "Computer Generation of Hypergeometric
    /// Random Variates" by VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
    ///  "COMPUTER GENERATION OF HYPERGEOMETRIC RANDOM VARIATES,"
    ///  JOURNAL OF STATISTICAL COMPUTATION AND SIMULATION,
    ///  22(1985), 2, 1985, 127-145.
    fn h2pec(kk: u64, nn1: u64, nn2: u64, rng: &mut ThreadRng) -> u64 {

        // input validation
        let (n1, n2) = match _validate_input(kk as f64, nn1 as f64, nn2 as f64) {
            Ok(v)=> v,
            Err(errstr) => panic!(errstr),
        };

        // variable setup
        let tn = n1 + n2;
        let k = if (kk + kk) as f64 >= tn {
            tn - kk as f64
        } else {
            kk as f64
        };
        let m = (k+1.0) * (n1+1.0) / (tn+2.0);

        let minjx:f64 = if k-n2 < 0.0 {
            0.0
        } else {
            k-n2
        };
        let maxjx:f64 = if n1 < k {
            n1
        } else {
            k
        };

        // generate random variate
        let ix: f64 = if minjx == maxjx {
            maxjx
        } else if m - minjx < 10.0 {
            // inverse transformation
            let w = if k < n2 {
                (CON + _afc(n2).unwrap() + _afc(n1+n2-k).unwrap()
                    - _afc(n2-k).unwrap() - _afc(n1+n2).unwrap()).exp()
            } else {
                (CON + _afc(n1).unwrap() + _afc(k).unwrap()
                    - _afc(k-n2).unwrap() - _afc(n1+n2).unwrap()).exp()
            };

            let mut lix;
            'label10: loop {
                let mut p = w;
                lix = minjx;
                let mut u = _rand(rng) * SCALE;
                'label20: loop {
                    if u>p {
                        u = u - p;
                        p = p* (n1 - lix) * (k-lix);
                        lix = lix + 1.0;
                        p = p / lix / (n2 - k + lix);
                        if lix > maxjx {
                            continue 'label10
                        }
                        continue
                    }
                    break 'label10
                }
            }
            lix
        } else {
            // h2pe
            let mut lix;

            let s = ((tn-k) * k * n1 * n2 / (tn -1.0)/tn/tn).sqrt();
            let d = (1.5 * s).trunc() + 0.5;
            let xl = (m - d + 0.5).trunc();
            let xr = (m + d + 0.5).trunc();
            let a = _afc(m).unwrap() + _afc(n1-m).unwrap() + _afc(k-m).unwrap()
                + _afc(n2-k+m).unwrap();
            let expon = a - _afc(xl).unwrap() - _afc(n1-xl).unwrap() 
                - _afc(k-xl).unwrap() - _afc(n2-k+xl).unwrap();
            let kl = expon.exp();
            let kr = (a - _afc(xr-1.0).unwrap() - _afc(n1-xr+1.0).unwrap()
                      - _afc(k-xr+1.0).unwrap() - _afc(n2-k+xr-1.0).unwrap()).exp();
            let lamdl = -1.0 * (xl * (n2 - k + xl)/(n1-xl+1.0)/(k-xl+1.0)).ln();
            let lamdr = -1.0 * ((n1-xr+1.0) * (k-xr+1.0)/xr/(n2-k+xr)).ln();

            let p1 = 2.0*d;
            let p2 = p1 + kl / lamdl;
            let p3 = p2 + kr / lamdr;


            'label30: loop {
                let mut reject = true;
                let mut u = _rand(rng) * p3;
                let mut v = _rand(rng);

                if u < p1 {
                    // rect region
                    lix = xl + u;
                } else if u <= p2 {
                    lix = xl + v.ln()/lamdl;
                    if lix < minjx {
                        continue;
                    }
                    v = v * (u-p1) * lamdl;
                } else {
                    lix = xr - v.ln()/lamdr;
                    if lix > maxjx {
                        continue;
                    }
                    v = v * (u-p2) * lamdr;
                }

                if m < 100.0 || lix <= 50.0 {
                    // explicit evaluation
                    let mut f = 1.0;
                    if m < lix {
                        let mut i = m+1.0;
                        while i < lix {
                            f = f * (n1 - i + 1.0) * (k - i + 1.0)
                                / (n1 - i) / (k-i);
                            i+=1.0;
                        }
                    } else if  m > lix  {
                        let mut i = lix+1.0;
                        while i < m {
                            f = f * i * (n2 - k + i) * (n1 - i) / (k-i);
                            i+=1.0;
                        }

                    }
                    if v <= f {
                        reject = false;
                    }
                } else {
                    // squeeze using upper and lower bounds
                    let y = lix;
                    let y1 = y+1.0;
                    let ym = y - m;
                    let yn = n1 - y + 1.0;
                    let yk = k - y + 1.0;
                    let nk = n2 - k + y1;
                    let r = -1.0 * ym / y1;
                    let s = ym / yn;
                    let t = ym / yk;
                    let e = -1.0 * t;
                    let g = yn * yk / (y1 * nk) - 1.0;
                    let dg = if g < 0.0 {
                        1.0 + g
                    } else {
                        1.0
                    };
                    let gu = g * (1.0 + g * (-0.5+g/3.0));
                    let gl = gu - 0.25 * ((g*g)*(g*g))/dg;
                    let xm = m + 0.5;
                    let xn = n1 - m + 0.5;
                    let xk = k - m + 0.5;
                    let nm = n2 - k + xm;
                    let ub = y * gu - m * gl + DELTAU +
                        xm * r * (1.0 + r * (-0.5+r/3.0)) +
                        xn * s * (1.0 + s * (-0.5+s/3.0)) +
                        xk * t * (1.0 + t * (-0.5+t/3.0)) +
                        nm * e * (1.0 + e * (-0.5+e/3.0));

                    let alv = v.ln();
                    if alv > ub {
                        // test upper bounds
                        reject = true;
                    } else {
                        // lower bounds
                        let mut dr = xm * ((r * r) * (r * r));
                        if r < 0.0 {
                            dr = dr / (1.0+r);
                        }
                        let mut ds = xn * ((s * s) * (s * s));
                        if s < 0.0 {
                            ds = ds / (1.0 + s);
                        }
                        let mut dt = xk * ((t * t) * (t * t));
                        if t < 0.0 {
                            dt = dt / (1.0 + t);
                        }
                        let mut de = nm * ((e * e) * (e * e));
                        if e < 0.0 {
                            de = de / (1.0 + e);
                        }
                        if alv < ub - 0.25 * (dr+ds+dt+de) + (y+m) * (gl-gu) - DELTAL {
                            reject = false;

                        } else {
                            if alv <= (a - _afc(lix).unwrap() - _afc(n1-lix).unwrap() - _afc(k-lix).unwrap() - _afc(n2 - k + lix).unwrap() ) {
                                reject = false;
                            } else {
                                reject = true;
                            }
                        }
                    }
                }
                if reject {
                    continue
                }
                break;
            }
            lix
        };

        // return variate
        if kk + kk >= tn.round() as u64 {
            if nn1 > nn2 {
                (kk as f64 - nn2 as f64 + ix).round() as u64
            } else {
                (nn1 as f64 - ix).round() as u64
            }
        } else {
            if nn1 > nn2 {
                (kk as f64 - ix) as u64
            } else {
                ix.round() as u64
            }
        }
    }


    #[test]
    fn _test_h2pec() {
        let mut rng = thread_rng();
        let _jx = h2pec(10, 100, 100, &mut rng);
    }


    fn _validate_input(kk:f64, nn1:f64, nn2:f64) -> Result<(f64,f64),String> {
        if nn1 < 0.0 || nn2 < 0.0 || kk < 0.0 || kk>(nn1+nn2) {
            Err(String::from("Invalid inputs"))
        } else {
            Ok(if nn1 >= nn2 {
                (nn2, nn1)
            } else {
                (nn1, nn2)
            })
        }
    }
    #[test]
    fn _test_validate_input() {
        assert!(
            match _validate_input(2.0, 1.0, 1.0) {
                Ok(_) => true,
                Err(_) => false,
            });
        assert!(
            match _validate_input(3.0, 1.0, 1.0) {
                Ok(_) => false,
                Err(_) => true,
            });
        assert!(
            match _validate_input(2.0, -1.0, 1.0) {
                Ok(_) => false,
                Err(_) => true,
            });
        assert!(
            match _validate_input(2.0, 1.0, -1.0) {
                Ok(_) => false,
                Err(_) => true,
            });
        assert!(
            match _validate_input(-2.0, 1.0, 1.0) {
                Ok(_) => false,
                Err(_) => true,
            });
    }


    // _afc -
    // Function to evaluate logarithm of the factorial i.
    // if i >= 7 use stirling's approximation
    fn _afc(i: f64) -> Option<f64> {
        if i < 0.0 {
            None
        } else {
            // if lte 7 use computed table
            match i.round() as u64 {
                0 => Some(0.0),
                1 => Some(0.0),
                2 => Some(0.6931471806),
                3 => Some(1.791759469),
                4 => Some(3.178053830),
                5 => Some(4.787491743),
                6 => Some(6.579251212),
                7 => Some(8.525161361),
                _ => Some(
                    (i+0.5) * i.ln() - i + 0.08333333333333/i
                    -0.00277777777777/i/i/i + 0.9189385332)
            }
        }
    }

    #[test]
    fn _test_afc() {
        assert_eq!(_afc(1.0).unwrap(), 0.0);
        assert_eq!(_afc(7.0).unwrap(), 8.525161361);
        assert_eq!(_afc(8.0).unwrap(), 10.604602878798048);
    }


    // rand - cryptographic random; uses range 0-1
    fn _rand(rng: &mut ThreadRng) -> f64 {
        let range = Uniform::new(0.0f64, 1.0f64);
        range.sample(rng)
    }

    #[test]
    fn _test_rand() {
        let mut rng = thread_rng();
        let result = _rand(&mut rng);
        assert!(result >=0.0 && result <=1.0);
    }
}



#[cfg(test)]
mod tests {

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
