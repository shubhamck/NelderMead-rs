# NelderMead-rs

Rust Implementtion of the popular Nelder Mead Algorithm for gradient-free optimization.

This implementation follows the explanation provided [here](https://codesachin.wordpress.com/2016/01/16/nelder-mead-optimization/).

Feel free to contribute and open issues.

## Usage

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

## To be Added
- Constrained Nelder Mead
- Parrallelized operations
- Status of Solution