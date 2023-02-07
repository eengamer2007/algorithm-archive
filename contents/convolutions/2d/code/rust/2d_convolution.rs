use std::cmp::max;

fn main() {
    
}

fn convolve_linear(signal: &[&[f64]], filter: &[&[f64]], outputsize: (isize, isize)) {
    let mut output: Vec<Vec<f64>> = Vec::with_capacity(outputsize.0 as usize);
    let mut sum: f64 = 0.0;

    for i in 0..outputsize.0 {
        for j in 0..outputsize.1 {
            for k in (max(0, i - filter.len() as isize) as usize)..=i {
                for l in (max(0, i - filter[0].len() as isize) as usize)..=i {
                    if k < 0 || k > filter.len() as isize || l < 0 || l > filter.len() as isize {
                        sum += 0.0;
                    } else {
                        sum += signal[k][l] * filter[i-k][j-l];
                    }
                }
            }
            output[i as usize][j as usize] = sum;
            sum = 0.0;
        }
    }

    
}

fn create_gaussian_kernel(kernel_size: usize) -> Vec<Vec<f64>> {
    let mut kernel: Vec<Vec<f64>> = Vec::with_capacity(kernel_size);

    // center must be offset by 0.5 to find the correct center
    let center: f64 = kernel_size as f64 * 0.5 + 0.5;

    let sigma = (0.1*kernel_size as f64).sqrt();

    for (i, x) in kernel.iter_mut().enumerate() {
        for (j, mut y) in x.iter_mut().enumerate() {
            y = &mut (
                -((i as f64-center).powf(2.0) + (j as f64-center).powf(2.0)) / ((2.0*sigma).powf(2.0))
            ).exp();
        }
    }

    normalize(kernel)
}


fn norm(array: &[Vec<f64>]) -> f64 {
    array.iter().map(|x| x.iter().sum::<f64>()).sum::<f64>().sqrt()
}

fn normalize(array: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let norm = norm(&array);

    let mut output: Vec<Vec<f64>> = Vec::with_capacity(array.len());

    for (x, j) in array.iter().enumerate() {
        for (y, i) in j.iter().enumerate() {
            output[x][y] = i / norm;
        }
    }

    output
}


fn create_sobel_operators() -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let sx: Vec<Vec<f64>> = vec![
        vec![1.0, 0.0, -1.0],
        vec![2.0, 0.0, -2.0],
        vec![1.0, 0.0, -1.0]
    ];

    let sx: Vec<Vec<f64>> = vec![
        vec![1.0, 2.0, 1.0],
        vec![0.0, 0.0, 0.0],
        vec![-1.0, -2.0, -1.0]
    ];

    (sx,sy)
}

fn sum_matrix_dimensions(array1: &[Vec<f64>], array2: &[Vec<f64>]) -> (isize, isize) {
    ((array1.len() + array2.len) as isize, (array1[0].len() + array2[0].len()) as isize)
}

fn compute_sobel(signal: &[Vec<f64>]) {
    (sx, sy) = create_sobel_operators();

    let gx = convolve_linear(signal, sx, sum_matrix_dimensions(signal, sx));
    let gy = convolve_linear(signal, sy, sum_matrix_dimensions(signal, sy));

    gx.iter().enumerate().map(|(x, i)| (x.powf(2.0) + gy[i].powf()).sqrt())
}
