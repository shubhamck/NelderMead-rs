extern crate nelder_mead_rs;

fn main() {
    // rosenbrock function (1 - x)^2 + 100 + (y - x^2 )^2
    let f = |v: &[f64]| -> f64 {
        let s = (1.0 - v[0]).powi(2) + 100.0 + (v[1] - v[0].powi(2)).powi(2);
        s
    };
    let initial_guess = vec![3.0, 4.0];
    // assert_eq!(nelder_mead_solve(&f, &initial_guess), [0.0, 0.0]);
    let solution = nelder_mead_rs::nm::nelder_mead_solve(&f, &initial_guess);
    println!("Solution is {:?}", solution);
}
