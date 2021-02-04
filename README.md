# NelderMead-rs

[![Build Status](https://api.travis-ci.com/shubhamck/NelderMead-rs.svg?branch=master)](https://travis-ci.com/github/shubhamck/NelderMead-rs)

Rust Implementtion of the popular Nelder Mead Algorithm for gradient-free optimization.

This implementation follows the explanation provided [here](https://codesachin.wordpress.com/2016/01/16/nelder-mead-optimization/).

Feel free to contribute and open issues.

## Usage

### Unconstrained Optimization

```rust
extern crate nelder_mead_rs;

fn main() {
    // Describe your cost function as a closure
    let f = |v: &[f64]| -> f64 {
        let mut s = 0.0;
        for x in v {
            s += x.powi(2);
        }
        s
    };
    // provide initial guess
    let initial_guess = vec![3.0, 4.0];
    // call solve
    let solution = nelder_mead_rs::nm::nelder_mead_solve(&f, &initial_guess);

}
```

### Constrained Optimization
```rust
extern crate nelder_mead_rs;

fn main() {
    // Describe your cost function as a closure
    let f = |v: &[f64]| -> f64 {
        let mut s = 0.0;
        for x in v {
            s += x.powi(2);
        }
        s
    };
    // provide initial guess
    let initial_guess = vec![3.0, 4.0];
    let lower_bounds = vec![1.0, 1.0];
    let upper_bounds = vec![100.0, 100.0];
    // call solve
    let solution = nelder_mead_rs::nm::nelder_mead_solve_constrained(&f, &initial_guess, &lower_bounds, &upper_bounds);

}
```

## To be Added
- Constrained Nelder Mead
- Parrallelized operations
- Status of Solution