fn gauleg(x1: f64, x2: f64, x: &mut [f64], w: &mut [f64]) {
    const EPS: f64 = 1.0e-14;
    let n = x.len();
    let m = (n + 1) / 2;
    let xm = 0.5 * (x2 + x1);
    let xl = 0.5 * (x2 - x1);
    for i in 0..m {
        let z = (3.141592654 * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        loop {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j as f64 + 1.0) * z * p2 - j as f64 * p3) / (j + 1) as f64;
            }
            let pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z -= p1 / pp;
            if (z - z1).abs() <= EPS { break; }
        }
        x[i] = xm - xl * z;
        x[n - 1 - i] = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
}

fn gaulag(x: &mut [f64], w: &mut [f64], alf: f64) {
    const MAXIT: i32 = 10;
    const EPS: f64 = 1.0e-14;
    let n = x.len();
    for i in 0..n {
        if i == 0 {
            let z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n as f64 + 1.8 * alf);
        } else if i == 1 {
            let z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n as f64);
        } else {
            let ai = i as f64 - 1.0;
            let z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf / (1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * alf);
        }
        for its in 0..MAXIT {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j + 1 + alf - z) * p2 - (j + alf) * p3) / (j + 1) as f64;
            }
            let pp = (n as f64 * p1 - (n + alf) as f64 * p2) / z;
            let z1 = z;
            z -= p1 / pp;
            if (z - z1).abs() <= EPS { break; }
        }
        if its >= MAXIT { panic!("too many iterations in gaulag"); }
        x[i] = z;
        w[i] = -((gammln(alf + n as f64) - gammln(n as f64)) / (pp * n as f64 * p2)).exp();
    }
}

fn gauher(x: &mut [f64], w: &mut [f64]) {
    const EPS: f64 = 1.0e-14;
    const PIM4: f64 = 0.7511255444649425;
    const MAXIT: i32 = 10;
    let n = x.len();
    let m = (n + 1) / 2;
    for i in 0..m {
        if i == 0 {
            let z = (2 * n + 1).sqrt() - 1.85575 * ((2 * n + 1) as f64).powf(-0.16667);
        } else if i == 1 {
            let z -= 1.14 * ((n as f64).powf(0.426) / z);
        } else if i == 2 {
            let z *= 1.86;
            z -= 0.86 * x[0];
        } else if i == 3 {
            let z *= 1.91;
            z -= 0.91 * x[1];
        } else {
            z = 2.0 * z - x[i - 2];
        }
        for its in 0..MAXIT {
            let mut p1 = PIM4;
            let mut p2 = 0.0;
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = z * (2.0 / (j + 1) as f64).sqrt() * p2 - ((j as f64) / (j + 1) as f64).sqrt() * p3;
            }
            let pp = ((2 * n) as f64).sqrt() * p2;
            let z1 = z;
            z -= p1 / pp;
            if (z - z1).abs() <= EPS { break; }
        }
        if its >= MAXIT { panic!("too many iterations in gauher"); }
        x[i] = z;
        x[n - 1 - i] = -z;
        w[i] = 2.0 / (pp * pp);
        w[n - 1 - i] = w[i];
    }
}

fn gaujac(x: &mut [f64], w: &mut [f64], alf: f64, bet: f64) {
    const MAXIT: i32 = 10;
    const EPS: f64 = 1.0e-14;
    let n = x.len();
    let mut an: f64 = 0.0;
    let mut bn: f64 = 0.0;
    let mut r1: f64 = 0.0;
    let mut r2: f64 = 0.0;
    let mut r3: f64 = 0.0;
    let mut z: f64 = 0.0;
    let mut temp: f64 = 0.0;
    let mut p1: f64 = 0.0;
    let mut p2: f64 = 0.0;
    let mut p3: f64 = 0.0;
    let mut pp: f64 = 0.0;
    let mut its: i32 = 0;
    let mut an: f64 = 0.0;
    let mut bn: f64 = 0.0;
    let alfbet: f64 = alf + bet;

    for i in 0..n {
        match i {
            0 => {
                an = alf / (n as f64);
                bn = bet / (n as f64);
                r1 = (1.0 + alf) * (2.78 / (4.0 + (n * n) as f64) + 0.768 * an / (n as f64));
                r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an * bn;
                z = 1.0 - r1 / r2;
            },
            1 => {
                r1 = (4.1 + alf) / ((1.0 + alf) * (1.0 + 0.156 * alf));
                r2 = 1.0 + 0.06 * ((n as f64) - 8.0) * (1.0 + 0.12 * alf) / (n as f64);
                r3 = 1.0 + 0.012 * bet * (1.0 + 0.25 * (alf).abs()) / (n as f64);
                z -= (1.0 - z) * r1 * r2 * r3;
            },
            2 => {
                r1 = (1.67 + 0.28 * alf) / (1.0 + 0.37 * alf);
                r2 = 1.0 + 0.22 * ((n as f64) - 8.0) / (n as f64);
                r3 = 1.0 + 8.0 * bet / ((6.28 + bet) * (n * n) as f64);
                z -= (x[0] - z) * r1 * r2 * r3;
            },
            n - 2 => {
                r1 = (1.0 + 0.235 * bet) / (0.766 + 0.119 * bet);
                r2 = 1.0 / (1.0 + 0.639 * ((n as f64) - 4.0) / (1.0 + 0.71 * ((n as f64) - 4.0)));
                r3 = 1.0 / (1.0 + 20.0 * alf / ((7.5 + alf) * (n * n) as f64));
                z += (z - x[n - 4]) * r1 * r2 * r3;
            },
            n - 1 => {
                r1 = (1.0 + 0.37 * bet) / (1.67 + 0.28 * bet);
                r2 = 1.0 / (1.0 + 0.22 * ((n as f64) - 8.0) / (n as f64));
                r3 = 1.0 / (1.0 + 8.0 * alf / ((6.28 + alf) * (n * n) as f64));
                z += (z - x[n - 3]) * r1 * r2 * r3;
            },
            _ => {
                z = 3.0 * x[i - 1] - 3.0 * x[i - 2] + x[i - 3];
            }
            
        }
    }
    alfbet = alf + bet;

    for its in 1..=MAXIT {
        temp = 2.0 + alfbet;
        p1 = (alf - bet + temp * z) / 2.0;
        p2 = 1.0;

        for j in 2..=n {
            p3 = p2;
            p2 = p1;
            temp = 2.0 * j as f64 + alfbet;
            let a = 2.0 * j as f64 * (j as f64 + alfbet) * (temp - 2.0);
            let b = (temp - 1.0) * (alf * alf - bet * bet + temp * (temp - 2.0) * z);
            let c = 2.0 * (j - 1 + alf) * (j - 1 + bet) * temp;
            p1 = (b * p2 - c * p3) / a;
        }

        let pp = (n as f64 * (alf - bet - temp * z) * p1 + 2.0 * ((n as f64) + alf) * ((n as f64) + bet) * p2) / (temp * (1.0 - z * z));
        let z1 = z;
        z = z1 - p1 / pp;

        if (z - z1).abs() <= EPS {
            break;
        }
    }

    if its > MAXIT {
        panic!("too many iterations in gaujac");
    }

    x[i] = z;

    let w[i] = (gammln(alf + n as f64) + gammln(bet + n as f64) - gammln((n + 1) as f64) - gammln(n + alfbet + 1.0)) * temp * 2.0_f64.powf(alfbet) / (pp * p2);
}