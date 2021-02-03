fn get_unit_vector(initial_guess: &[f64]) -> Vec<f64> {
    // get the norm
    let sq_norm: f64 = initial_guess.iter().map(|x| x * x).sum();
    let norm = sq_norm.sqrt();

    // divide each element by norm
    let unit_vector: Vec<f64> = initial_guess.iter().map(|x| x / norm).collect();
    unit_vector
}

fn get_initial_simplex(initial_guess: &[f64]) -> Vec<Vec<f64>> {
    let dim = initial_guess.len();
    let unit_vector = get_unit_vector(initial_guess);
    let mut simplex: Vec<Vec<f64>> = Vec::new();
    let initial_guess_vec = initial_guess.to_vec();

    const SCALE: f64 = 0.05;

    simplex.push(initial_guess_vec.clone());

    for i in 0..dim {
        let mut v: Vec<f64> = initial_guess_vec.clone();
        v[i] += unit_vector[i] * SCALE;
        simplex.push(v);
    }
    simplex
}

fn calculate_costs(simplex: &[Vec<f64>], cost_function: &dyn Fn(&[f64]) -> f64) -> Vec<f64> {
    let costs: Vec<f64> = simplex.iter().map(|x| cost_function(x)).collect();
    return costs;
}

fn calculate_ordering(costs: &[f64]) -> Vec<usize> {
    let samples = costs.len();
    let mut ordering: Vec<usize> = (0..samples).collect();
    ordering.sort_by(|i, j| costs[*i].partial_cmp(&costs[*j]).unwrap());
    return ordering;
}

fn reorder(simplex: &mut Vec<Vec<f64>>, costs: &mut Vec<f64>) {
    let ordering = calculate_ordering(costs);

    let mut reordered_cost: Vec<f64> = Vec::new();
    let mut reordered_simplex: Vec<Vec<f64>> = Vec::new();

    for i in ordering {
        reordered_cost.push(costs[i]);
        reordered_simplex.push(simplex[i].clone());
    }

    *costs = reordered_cost;
    *simplex = reordered_simplex;
}

fn calculate_centroid(sorted_simplex: &[Vec<f64>]) -> Vec<f64> {
    let dim = sorted_simplex[0].len();
    let mut centroid = vec![0.0; dim];
    for i in 0..dim {
        let mut c = 0.0;
        for j in 0..sorted_simplex.len() {
            c += sorted_simplex[j][i];
        }
        c = c / (sorted_simplex.len()) as f64;
        centroid[i] = c;
    }
    centroid
}

fn calculate_reflected_point(pivot: &[f64], point: &[f64], reflection_factor: f64) -> Vec<f64> {
    let mut reflected_point = vec![0.0; pivot.len()];
    for (i, (aval, bval)) in pivot.iter().zip(point).enumerate() {
        reflected_point[i] = aval + reflection_factor * (aval - bval);
    }
    reflected_point
}

fn calculate_expanded_point(pivot: &[f64], point: &[f64], expansion_factor: f64) -> Vec<f64> {
    let mut expanded_point = vec![0.0; pivot.len()];
    for (i, (aval, bval)) in pivot.iter().zip(point).enumerate() {
        expanded_point[i] = aval + expansion_factor * (bval - aval);
    }
    expanded_point
}

fn calculate_contracted_point(pivot: &[f64], point: &[f64], contraction_factor: f64) -> Vec<f64> {
    let mut contracted_point = vec![0.0; pivot.len()];
    for (i, (aval, bval)) in pivot.iter().zip(point).enumerate() {
        contracted_point[i] = aval + contraction_factor * (bval - aval);
    }
    contracted_point
}

fn shrinkage_contraction(sorted_simplex: &mut Vec<Vec<f64>>, shrinkage_contraction_factor: f64) {
    // pivot is the best point in the sorted simplex
    let pivot = sorted_simplex[0].clone();
    let samples = sorted_simplex.len();
    for i in 1..samples {
        let contracted_point =
            calculate_contracted_point(&pivot, &sorted_simplex[i], shrinkage_contraction_factor);

        sorted_simplex[i] = contracted_point;
    }
}

pub fn nelder_mead_solve(cost_function: &dyn Fn(&[f64]) -> f64, initial_guess: &[f64]) -> Vec<f64> {
    let mut initial_simplex = get_initial_simplex(initial_guess);
    let num_samples = initial_simplex.len();
    let best_index: usize = 0;
    let worst_index: usize = num_samples - 1;
    let second_worst_index: usize = num_samples - 2;
    let mut iter = 0;
    const MAX_ITER: usize = 5000;
    const REFLECTION_FACTOR: f64 = 1.0;
    const EXPANSION_FACTOR: f64 = 2.0;
    const CONTRACTION_FACTOR: f64 = 0.5;
    const SHRINKAGE_CONTRACTION_FACTOR: f64 = 0.5;
    const BEST_COST_THRESHOLD: f64 = 1.0e-9;

    loop {
        let mut costs = calculate_costs(&initial_simplex, cost_function);

        // get centroid
        reorder(&mut initial_simplex, &mut costs);
        let best_cost = costs[best_index];
        let worst_cost = costs[worst_index];
        let second_worst_cost = costs[second_worst_index];
        if best_cost.abs() < BEST_COST_THRESHOLD {
            return initial_simplex[0].clone();
        }

        let centroid = calculate_centroid(&initial_simplex[..num_samples - 1]);

        // get reflected worst point
        let reflected_point =
            calculate_reflected_point(&centroid, &initial_simplex[worst_index], REFLECTION_FACTOR);

        let reflected_point_cost = cost_function(&reflected_point);

        // Three cases
        // 1. reflected cost is better than the best cost -> calculate_expanded_point
        // 2. reflected cost is better than second worst cost -> replace worst point wth reflected
        //    point
        // 3. reflected cost is worst than current worst
        if reflected_point_cost < best_cost {
            // reflected point is better than best point
            let expanded_point =
                calculate_expanded_point(&centroid, &reflected_point, EXPANSION_FACTOR);
            // replace worst point with better of reflected_point or expanded point
            let expanded_point_cost = cost_function(&expanded_point);

            if expanded_point_cost < reflected_point_cost {
                initial_simplex[worst_index] = expanded_point;
            } else {
                initial_simplex[worst_index] = reflected_point;
            }
        } else if reflected_point_cost < second_worst_cost {
            // reflected point is worst than best but better than second worst
            initial_simplex[worst_index] = reflected_point;
        } else if reflected_point_cost < worst_cost {
            // reflected point is worse than second worst but better than worst
            let contracted_point = calculate_contracted_point(
                &centroid,
                &initial_simplex[worst_index],
                CONTRACTION_FACTOR,
            );

            let contracted_point_cost = cost_function(&contracted_point);
            if contracted_point_cost < worst_cost {
                // contracted point is better than worst
                initial_simplex[worst_index] = contracted_point;
            } else {
                shrinkage_contraction(&mut initial_simplex, SHRINKAGE_CONTRACTION_FACTOR);
            }
        } else {
            shrinkage_contraction(&mut initial_simplex, SHRINKAGE_CONTRACTION_FACTOR);
        }

        iter += 1;

        if iter > MAX_ITER {
            break;
        }
    }
    return initial_simplex[0].clone();
}

#[cfg(test)]
mod tests {

    use super::calculate_centroid;
    use super::calculate_contracted_point;
    use super::calculate_costs;
    use super::calculate_expanded_point;
    use super::calculate_ordering;
    use super::calculate_reflected_point;
    use super::get_initial_simplex;
    use super::get_unit_vector;
    use super::nelder_mead_solve;
    use super::reorder;
    use super::shrinkage_contraction;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_solve() {
        let f = |v: &[f64]| -> f64 {
            let mut s = 0.0;
            for x in v {
                s += x.powi(2);
            }
            s
        };
        let initial_guess = vec![3.0, 4.0];
        // assert_eq!(nelder_mead_solve(&f, &initial_guess), [0.0, 0.0]);
        let solution = nelder_mead_solve(&f, &initial_guess);
        let expected_solution = vec![0.0, 0.0];
        for i in 0..initial_guess.len() {
            assert_approx_eq!(solution[i], expected_solution[i], 0.001);
        }
    }

    #[test]
    fn test_unit_vector() {
        let vector = vec![3.0, 4.0];
        assert_eq!(get_unit_vector(&vector), [0.6, 0.8]);
    }

    #[test]
    fn test_initial_simplex() {
        let vector = vec![3.0, 4.0];
        assert_eq!(
            get_initial_simplex(&vector),
            [[3.0, 4.0], [3.0 + 0.05 * 0.6, 4.0], [3.0, 4.0 + 0.05 * 0.8]]
        );
    }

    #[test]
    fn test_calculate_cost() {
        let f = |v: &[f64]| -> f64 {
            let mut s = 0.0;
            for x in v {
                s += x;
            }
            s
        };
        let simplex = vec![vec![1.0, 3.0], vec![4.0, 5.0]];
        let costs = calculate_costs(&simplex, &f);
        assert_eq!(costs, [4.0, 9.0]);
    }

    #[test]
    fn test_calculate_centroid() {
        let sorted_simplex = vec![vec![3.0, 4.0], vec![5.0, 6.0], vec![1.0, 2.0]];
        assert_eq!(calculate_centroid(&sorted_simplex), [3.0, 4.0]);
    }

    #[test]
    fn test_calculate_reflected_point() {
        let pivot = vec![0.0, 1.0];
        let point = vec![1.0, 1.0];
        assert_eq!(calculate_reflected_point(&pivot, &point, 1.0), [-1.0, 1.0]);
    }

    #[test]
    fn test_calculate_expanded_point() {
        let pivot = vec![0.0, 1.0];
        let point = vec![1.0, 1.0];
        assert_eq!(calculate_expanded_point(&pivot, &point, 2.0), [2.0, 1.0]);
    }

    #[test]
    fn test_calculate_contracted_point() {
        let pivot = vec![0.0, 1.0];
        let point = vec![1.0, 1.0];
        assert_eq!(calculate_contracted_point(&pivot, &point, 0.5), [0.5, 1.0]);
    }

    #[test]
    fn test_calculate_ordering() {
        let costs = vec![4.0, 2.0, 1.0, 3.0];
        assert_eq!(calculate_ordering(&costs), [2, 1, 3, 0]);
    }

    #[test]
    fn test_reorder() {
        let mut simplex = vec![vec![3.0, 4.0], vec![5.0, 6.0], vec![1.0, 2.0]];
        let mut costs = vec![3.0, 1.0, 2.0];
        reorder(&mut simplex, &mut costs);
        assert_eq!(costs, [1.0, 2.0, 3.0]);
        assert_eq!(
            simplex,
            vec![vec![5.0, 6.0], vec![1.0, 2.0], vec![3.0, 4.0]]
        );
    }

    #[test]
    fn test_shrinkage_contraction() {
        let mut simplex = vec![vec![3.0, 4.0], vec![5.0, 6.0], vec![1.0, 2.0]];
        shrinkage_contraction(&mut simplex, 0.5);
        assert_eq!(
            simplex,
            vec![vec![3.0, 4.0], vec![4.0, 5.0], vec![2.0, 3.0]]
        )
    }

    #[test]
    fn test_rosenbrock() {
        let f = |v: &[f64]| -> f64 {
            let s = (1.0 - v[0]).powi(2) + 100.0 + (v[1] - v[0].powi(2)).powi(2);
            s
        };
        let initial_guess = vec![3.0, 4.0];
        // assert_eq!(nelder_mead_solve(&f, &initial_guess), [0.0, 0.0]);
        let solution = nelder_mead_solve(&f, &initial_guess);
        let expected_solution = vec![1.0, 1.0];
        for i in 0..initial_guess.len() {
            assert_approx_eq!(solution[i], expected_solution[i], 0.001);
        }
    }
}
