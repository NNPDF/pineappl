#include <math.h>
#include <float.h>

double fx2(double y) {
    double x, yp, deltap, delta, deriv;
    int i;

    yp = y;
    deltap = DBL_MAX;

    for (i = 0; i <= 10; i=i+1) {
        x = exp(-yp);
        delta = (1.0 - x) * (-5.0) + (y - yp);

        if (fabs(delta) < 1e-15 && delta >= deltap) {
            return x;
        }

        deriv = -5.0 * x - 1.0;
        yp = yp - delta / deriv;
        deltap = delta;
    }

    return -1;
}
