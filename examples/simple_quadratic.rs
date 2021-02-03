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
    // assert_eq!(nelder_mead_solve(&f, &initial_guess), [0.0, 0.0]);
    let solution = nelder_mead_rs::nm::nelder_mead_solve(&f, &initial_guess);

    println!("Solution is {:?}", solution);
}
