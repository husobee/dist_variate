extern crate crypto;
// randvariate - module that provides functions to generate random variates
pub mod hgd {

    use crypto::symmetriccipher::SynchronousStreamCipher;

    pub struct Prng {
        pub cipher: Box<SynchronousStreamCipher + 'static>,
    }

    // coins is a bit string...
    // make sure the bit string is length 32
    // start with out = 0

    impl Prng {
        pub fn new(cipher: Box<SynchronousStreamCipher + 'static>) -> Prng {
            // new takes a derived key, and performs the encryption
            Prng{
                cipher: cipher,
            }
        }

        pub fn draw(&mut self) -> f64 {
            let mut coins = vec![0;16];
            self.cipher.process(&[0;16], &mut coins);

            let s: u32 = (coins[0] as u32) << 24 |
                (coins[1] as u32) << 16 | (coins[2] as u32) << 8 |
                (coins[3] as u32);

            let ret = (s as f64)/(u32::max_value() as f64);
            ret
        }
    }


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
    pub fn h2pec(kk: u64, nn1: u64, nn2: u64, rng: &mut Prng) -> u64 {

        // input validation
        let (n1, n2) = match _validate_input(kk as f64, nn1 as f64, nn2 as f64) {
            Ok(v)=> v,
            Err(errstr) => panic!(errstr),
        };

        // variable setup
        let tn = n1 + n2;
        let k = if (kk as f64 + kk as f64) >= tn {
            tn - kk as f64
        } else {
            kk as f64
        };
        let m = ((k+1.0) * (n1+1.0) / (tn+2.0)).trunc();

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
                (CON + _afc(n2) + _afc(n1+n2-k)
                    - _afc(n2-k) - _afc(n1+n2)).exp()
            } else {
                (CON + _afc(n1) + _afc(k)
                    - _afc(k-n2) - _afc(n1+n2)).exp()
            };

            let mut lix;
            'label10: loop {
                let mut p = w;
                lix = minjx;
                let mut u = rng.draw() * SCALE;
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

            let a = _afc(m) + _afc(n1-m) + _afc(k-m)
                + _afc(n2-k+m);


            let kl = (
                a - _afc(xl)
                - _afc(n1-xl)
                - _afc(k-xl) 
                - _afc(n2-k+xl)).exp();

            let kr = (a - _afc((xr-1.0).trunc()) - _afc((n1-xr+1.0).trunc())
                      - _afc((k-xr+1.0).trunc()) - _afc((n2-k+xr-1.0).trunc())).exp();
            let lamdl = -1.0 * (xl * (n2 - k + xl)/(n1-xl+1.0)/(k-xl+1.0)).ln();
            let lamdr = -1.0 * ((n1-xr+1.0) * (k-xr+1.0)/xr/(n2-k+xr)).ln();

            let p1 = d+d;
            let p2 = p1 + kl / lamdl;
            let p3 = p2 + kr / lamdr;


            let mut reject = true;

            'label30: loop {
                let mut u = rng.draw() * p3;
                let mut v = rng.draw();

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
                                / (n2 - k + i) / i;
                            i+=1.0;
                        }
                    } else if  m > lix  {
                        let mut i = lix+1.0;
                        while i < m {
                            f = f * i * (n2 - k + i) / (n1 - i) / (k-i);
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
                    let e = -1.0 * ym / nk;
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
                            if alv <= (a - _afc(lix) - _afc(n1-lix) - _afc(k-lix) - _afc(n2 - k + lix) ) {
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
        let seed: Vec<u8> = vec![0, 0, 0, 42];
        let mut rng = Prng::new(&seed);
        let M = u32::max_value() as u64;
        let N = u64::max_value();
        let y = 0.0 + ((M as f64)/2.0).ceil();
        let jx = h2pec(y as u64, M, N-M, &mut rng);
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
    fn _afc(i: f64) -> f64 {
        // if lte 7 use computed table
        match i.round() as u64 {
            0 => 0.0,
            1 => 0.0,
            2 => 0.6931471806,
            3 => 1.791759469,
            4 => 3.178053830,
            5 => 4.787491743,
            6 => 6.579251212,
            7 => 8.525161361,
            _ => (i+0.5) * i.ln() - i + 0.08333333333333/i
                -0.00277777777777/i/i/i + 0.9189385332
        }
    }

    #[test]
    fn _test_afc() {
        assert_eq!(_afc(1.0), 0.0);
        assert_eq!(_afc(7.0), 8.525161361);
        assert_eq!(_afc(8.0), 10.604602878798048);
    }




    pub fn rhyper(kk:u64, nn1:u64, nn2:u64, prng: &mut Prng) -> u64 {
        if kk > 10 {
            hrua(kk, nn1, nn2, prng)
        } else {
            hyp(kk, nn1, nn2, prng)
        }
    }

    use std::cmp;

    fn hyp(kk:u64, nn1:u64, nn2:u64, prng: &mut Prng) -> u64 {
        let d1 = nn2 + nn1 - kk;
        let d2 = cmp::min(nn2, nn1);

        let mut Y = d2;
        let mut K = kk;

        while Y>0 {
            let U = prng.draw();
            Y = Y - (U + Y as f64 /(d1+K) as f64).floor() as u64;
            K = K - 1;
            if K == 0 {
                break;
            }
        }

        let Z = (d2 - Y) as u64;

        if nn1 > nn2 {
            kk - Z
        } else {
            Z
        }
    }

    fn hrua(kk:u64, nn1:u64, nn2:u64, prng: &mut Prng) -> u64 {
        const D1:f64 = 1.7155277699214135;
        const D2:f64 = 0.8989161620588988;

        let mut Z: u64 = 0;

        let mingoodbad = cmp::min(nn1, nn2);
        let maxgoodbad = cmp::max(nn1, nn2);
        let popsize = nn1+nn2;
        let m = cmp::min(kk, popsize - kk);

        let d4 = (mingoodbad/popsize) as f64;
        let d5 = 1.0 - d4;
        let d6 = m as f64 * d4 + 0.5;
        let d7 = ( (popsize - m) as f64 * kk as f64 * d4 * d5 / (popsize -1) as f64 +0.5).sqrt();
        let d8 = D1 * d7 + D2;
        let d9 = ((m+1)*(mingoodbad+1)/(popsize +2));
        let d10 = loggam(d9+1) + loggam(mingoodbad-d9+1) + loggam(m-d9+1) + loggam(maxgoodbad -m +d9+1);
        let d11 = cmp::min(cmp::min(m, mingoodbad), (d6+16.0*d7).floor() as u64);


        loop {
            let X = prng.draw();
            let Y = prng.draw();
            let W = d6 + d8 * (Y -0.5)/X;

            if W < 0.0 || W >= d11 as f64 {
                continue;
            }

            let Z = W.floor() as u64;
            let T = d10 - (loggam(Z+1) + loggam(mingoodbad - Z + 1) + loggam(m -Z+1)
                + loggam(maxgoodbad - m + Z +1));

            if (X * (4.0-X) -3.0) <= T {
                break;
            }

            if X * (X-T) >= 1.0 {
                continue;
            }

            if 2.0 * X.ln() <= T {
                break;
            }
        }

        if nn1 > nn2 {
            Z = m - Z;
        }

        if m < kk {
            Z = nn1 -Z
        }

        Z
    }
    use std::f64::consts;

    fn loggam(x:u64) -> f64 {

        let a: Vec<f64> = vec![
            8.333333333333333e-02, -2.777777777777778e-03,
            7.936507936507937e-04, -5.952380952380952e-04,
            8.417508417508418e-04, -1.917526917526918e-03,
            6.410256410256410e-03, -2.955065359477124e-02,
            1.796443723688307e-01, -1.39243221690590e+00];


        let mut x0 = x as f64;
        let mut n:u64 = 0;

        if x == 1 || x == 2 {
            return 0.0;
        } else if x <= 7 {
            let n = (7 - x);
            x0 = (x + n) as f64;
        }

        let x2 = 1.0 / (x0*x0);
        let xp = 2.0 * consts::PI;
        let mut gl0 = a[9];


        for v in (0..8).rev() {
            gl0 = gl0 * x2;
            gl0 = gl0 + a[v];
        }

        let mut gl = gl0/x0 as f64 + 0.5 * (xp as f64).ln() + (x0-0.5) * x0.ln() - x0;

        if x <= 7 {
            for i in (1..n+1) {
                gl = gl - (x0 as f64 - 1.0).ln();
                x0 = x0 - 1.0;
            }
        }
        gl
    }























    // rand - cryptographic random; uses range 0-1
    fn _rand(rng: &mut Prng) -> f64 {
        rng.draw()
    }

    #[test]
    fn _test_rand() {
        /*
        let seed: Vec<u8> = vec![255, 255, 255, 255];
        let mut rng = Prng::new(&seed);
        let result = _rand(&mut rng);
        assert!(result == 1.0);

        let seed: Vec<u8> = vec![0, 0, 0, 0];
        let mut rng = Prng::new(&seed);
        let result = _rand(&rng);
        assert!(result == 0.0);

        let seed: Vec<u8> = vec![0, 0, 0, 42];
        let rng = Prng::new(&seed);
        let result = _rand(&rng);
        assert!(result >=0.0 && result <=1.0);
        */
    }
}



#[cfg(test)]
mod tests {

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
