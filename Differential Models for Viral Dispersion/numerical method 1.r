Green.Vanishes = function(L, x, y, xi, eta) {
    return((x == 0 | x == L |
        xi == 0 | xi == L |
        (y > eta & (eta == 0 | y == L)) | (y <= eta & (y == 0 | eta == L))));
};

Green = function(L, x, y, xi, eta, precision = 1e-3) {
    result = 0;
    n = 1;
    if (Green.Vanishes(L, x, y, xi, eta)) return(0);
    while (T) {
        H = ifelse(
            y > eta,
            sinh(n * pi * eta / L) / (n * pi * sinh(n * pi)) * sinh(n * pi * (L - y) / L),
            sinh(n * pi * y / L) / (n * pi * sinh(n * pi)) * sinh(n * pi * (L - eta) / L)
        );
        summand = H * sin(n * pi * x / L) * sin(n * pi * xi / L);
        if (result != 0 & abs(summand / result) < precision) {
            return(result * 2);
        } else {
            result = result + summand;
            n = n + 1;
        }
    }
};

c.Solve = function(x, y, alpha.func, L, precision = 5 * 1e-4) {
    if (x == 0 | x == L | y == 0 | y == L) return(0);
    n = 4;
    last.integral = 0;
    while(T) {
        n.partition = 2^n;
        dxi = L / n.partition;
        deta = L / n.partition;
        xi.seq = seq(from = 0, to = L - dxi, by = dxi);
        eta.seq = seq(from = 0, to = L - deta, by = deta);
        dA = dxi * deta;
        integral = 0;
        for (xi in xi.seq) {
            for (eta in eta.seq) {
                alpha = alpha.func(xi + dxi / 2, eta + deta / 2);
                if (alpha != 0) {
                    integral = integral + alpha *
                        Green(L, x, y, xi + dxi / 2, eta + deta / 2) * dA;
                }
            }
        }
        if (abs(integral - last.integral) < precision) {
        #if (last.integral != 0 & abs((integral - last.integral) / last.integral) < precision) {
            return(integral);
        } else {
            last.integral = integral;
            n = n + 1;
        }
    }
};

# Accepts a sim.object which is a list of the form:
### sim.object$u: an array of viral count in each cell.
### sim.object$alpha: an array of viral creation rate in each cell.
### sim.object$kappa: a real value for diffusivity.
Iterate.1D.No.Wall = function(sim.object.1D, precision = 1e-8) {
    u = sim.object.1D$u;
    alpha = sim.object.1D$alpha;
    kappa = sim.object.1D$kappa;
    n.interval = length(u);
    last.u = rep(0, n.interval);
    while (T) {
        du = rep(0, n.interval);
        du[1] = alpha[1] - 2 * kappa * (u[1] - u[2] / 2);
        du[n.interval] = alpha[n.interval] - 2 * kappa * (u[n.interval] - u[n.interval - 1] / 2);
        du[2:(n.interval - 1)] = alpha[2:(n.interval - 1)] -
            2 * kappa * (u[2:(n.interval - 1)] - (u[1:(n.interval - 2)] + u[3:n.interval]) / 2);
        if (sum(du == 0) == 0 & max(abs(du / u)) < precision) {
            return(list(
                u = u,
                alpha = alpha,
                kappa = kappa
            ));
        } else {
            last.u = u;
            u = u + du;
        };
    }
};
    
Iterate.1D.Left.Wall = function(sim.object, precision = 1e-8) {
    u = sim.object$u;
    alpha = sim.object$alpha;
    kappa = sim.object$kappa;
    n.interval = length(u);
    last.u = rep(0, n.interval);
    while (T) {
        du = rep(0, n.interval);
        du[1] = alpha[1] - kappa * (u[1] - u[2]);
        du[n.interval] = alpha[n.interval] - 2 * kappa * (u[n.interval] - u[n.interval - 1] / 2);
        du[2:(n.interval - 1)] = alpha[2:(n.interval - 1)] -
            2 * kappa * (u[2:(n.interval - 1)] - (u[1:(n.interval - 2)] + u[3:n.interval]) / 2);
        if (sum(du == 0) == 0 & max(abs((du) / u)) < precision) {
            return(list(
                u = u,
                alpha = alpha,
                kappa = kappa
            ));
        } else {
            last.u = u;
            u = u + du;
        };
    }
};



Generate.c.Field = function(L, alpha.func, n.partition = 30) {
    mesh = L / (n.partition - 1);
    x.seq = seq(from = 0, to = L, by = mesh);
    y.seq = seq(from = 0, to = L, by = mesh);
    c.field = matrix(nrow = length(y.seq), ncol = length(x.seq));
    for (r in 1:length(y.seq)) {
        for (c in 1:length(x.seq)) {
            c.field[r, c] = c.Solve(x.seq[c], y.seq[r], alpha.func, L);
            print(sprintf("(%d, %d) done.", r, c));
        }
    }
    return(c.field);
}

# 1-dimensional scenarios.

# 1-D no wall, constant source.
sim.1 = Iterate.1D.No.Wall(list(
    u = rep(0, 100),
    alpha = rep(1, 100),
    kappa = 0.1
), precision = 1e-10);
plot(sim.1$u, col = c("red"), cex = 0.5);
x1.theoretical = 0:99 + 0.5;
c1.theoretical = - 1 / 0.2 * (x1.theoretical^2 - 100 * x1.theoretical);
lines(c1.theoretical);

sim.2 = Iterate.1D.Left.Wall(list(
    u = rep(0, 100),
    alpha = rep(1, 100),
    kappa = 0.1
), precision = 1e-10);
plot(sim.2$u, col = c("red"), cex = 0.5);
x2.theoretical = 0:99 + 0.5;
c2.theoretical = 1 / 0.2 * (100^2 - x2.theoretical^2);
lines(c2.theoretical);

# 1-D no wall, piecewise source.
sim.3 = Iterate.1D.No.Wall(list(
    u = rep(0, 100),
    alpha = c(rep(0, 20), rep(1, 20), rep(0, 20), rep(1, 20), rep(0, 20)),
    kappa = 0.1
), precision = 1e-10);
plot(sim.3$u, col = c("red"), cex = 0.5);
x3.sec1 = 0:19 + 0.5;
x3.sec2 = 20:39 + 0.5;
x3.sec3 = 40:59 + 0.5;
x3.sec4 = 60:79 + 0.5;
x3.sec5 = 80:99 + 0.5;
c3.theoretical = c(
    200 * x3.sec1,
    -10 * (0.5 * (x3.sec2^2) - 40 * x3.sec2 + 200),
    0 * x3.sec3 + 6000,
    -10 * (0.5 * (x3.sec4^2) - 60 * x3.sec4 + 1200),
    10 * (2000 - 20 * x3.sec5)
);
lines(c3.theoretical);
    
library(plot3D);

Persp.c.Field = function(c.field) {
    persp3D(1:ncol(c.field), 1:nrow(c.field), c.field, xlab = "x", ylab = "y", zlab = "c(x,y)");
};

c1.field = matrix(data = NA, nrow = 51, ncol = 51);
R = 25;
for (r in 1:nrow(c1.field)) {
    for (c in 1:ncol(c1.field)) {
        x = c - 26; y = 26 - r;
        if (x^2 + y^2 > R^2) {
            c1.field[r, c] = 0;
        } else {
            c1.field[r, c] = -(x^2 + y^2) + R^2;
        }
    }
}
c1.field;
Persp.c.Field(c1.field);

        

const.c.field = Generate.c.Field(L = 1, function(x, y) 1, n.partition = 50);
Persp.c.Field(const.c.field);

regional.alpha = function(x, y) {
    if ((x > 1/3 & x < 2/3) | (y > 1/3 & y < 2/3)) {
        return(0);
    } else {
        return(1);
    };
};
regional.c.field = Generate.c.Field(L = 1, regional.alpha, n.partition = 50);
Persp.c.Field(regional.c.field);

