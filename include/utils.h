#ifndef UTILS_H
#define UTILS_H

// Linear interpolation
double lerp(double a, double b, double t) {
    return a + (b-a)*t;
}

// Bilinear interpolation
//  Return: lerp(lerp(c00, c10, tx), lerp(c01, c11, tx), ty)
double blerp(double c00, double c10, double c01, double c11, double tx, double ty) {
    double a = c00 + (c10-c00)*tx;
    double b = c01 + (c11-c01)*tx;
    return (a + (b-a)*ty);
}

#endif