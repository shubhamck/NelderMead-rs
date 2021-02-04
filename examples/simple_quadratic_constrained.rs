extern crate nelder_mead_rs;

fn main() {
    // function is x^2 + y^2
    let f = |v: &[f64]| -> f64 {
        let mut s = 0.0;
        for x in v {
            s += x.powi(2);
        }
        s
    };
    let initial_guess = vec![3.0, 4.0];
    let lower_bounds = vec![1.0, 1.0];
    let upper_bounds = vec![100.0, 100.0];
    let solution = nelder_mead_rs::nm::nelder_mead_solve_constrained(
        &f,
        &initial_guess,
        &lower_bounds,
        &upper_bounds,
    );

    println!("Solution is {:?}", solution);
}
