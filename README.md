# randvariate

This code implements the H2PEC algorithm in rust.  It is based on the
following work:

```
VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER
"COMPUTER GENERATION OF HYPERGEOMETRIC RANDOM VARIATES,"
JOURNAL OF STATISTICAL COMPUTATION AND SIMULATION,
22(1985), 2, 1985, 127-145.
```

## usage

```rust
extern crate randvariate;
extern crate rand;

use randvariate::h2pec;


let mut rng = rand::thread_rng();
let jx = h2pec(100, 10000, 20000, &rng);

// jx will be the number of white balls drawn
```
