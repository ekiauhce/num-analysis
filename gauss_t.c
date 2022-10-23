#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "gauss.c"

struct method {
    void (*func)(struct matrix*);
    char* name;
};

void run_tests(struct method* method);
void _round_test();


const int VARIANT_3 = 3;
const int VARIANT_16 = 16;

int main(int argc, char const *argv[]) {
    mkdir("results", 0777);
    _round_test();

    struct method methods[] = {
        {calculate,               "calculate"},
        {calculate_with_pivoting, "calculate_with_pivoting"}
    };

    for (int i = 0; i < 2; i++) {
        run_tests(&methods[i]);
    }

    return 0;
}

double _round(double value, int places) {
    if (places == -1) {
        return value;
    }
    return round(value * pow(10, places)) / pow(10, places);
}

double** get_matrix_with_sqrts(int places, int variant) {
    int n = 4;
    double** a = malloc(sizeof(double*) * n);

    for (int i = 0; i < n; i++) {
        a[i] = malloc(sizeof(double) * n);
    }

    // https://miit.ru/content/Содержимое.pdf?id_vf=13922
    if (variant == VARIANT_3) {
        a[0][0] = 1.;                               a[0][1] = 1/_round(sqrt(2.), places);       a[0][2] =  1/_round(sqrt(2.), places); a[0][3] =  0;
        a[1][0] = 7./12;                            a[1][1] = 4./5;                             a[1][2] = -1./3;                       a[1][3] =  1./2;
        a[2][0] = 2./3;                             a[2][1] = _round(sqrt(2), places)/3;        a[2][2] = -1./3;                       a[2][3] = -1./_round(sqrt(3), places);
        a[3][0] = (1 + _round(sqrt(2), places))/2.; a[3][1] = (2 + _round(sqrt(2), places))/4.; a[3][2] =  0;                          a[3][3] = -_round(sqrt(3), places)/(2*_round(sqrt(2), places));
    } else if (variant == VARIANT_16) {
        a[0][0] = _round(sqrt(5), places);  a[0][1] =  2*_round(sqrt(3), places); a[0][2] =  0.5; a[0][3] = -1.4;
        a[1][0] = -_round(sqrt(3), places); a[1][1] = -6/_round(sqrt(5), places); a[1][2] =  3.2; a[1][3] =  2.5;
        a[2][0] = -1.1;                     a[2][1] =  4.2;                       a[2][2] = -3.2; a[2][3] =  0.8;
        a[3][0] =  7.9;                     a[3][1] = -5.2;                       a[3][2] =  0.3; a[3][3] =  2.9;
    }
    return a;
}

double* get_b_for_sqrts_case(int variant) {
    int n = 4;
    double* b = malloc(sizeof(double) * n);

    if (variant == VARIANT_3) {
        b[0] = 0.15;
        b[1] = 0.09;
        b[2] = -0.05;
        b[3] = 0.27;
    } else if (variant == VARIANT_3) {
        b[0] = 4.5;
        b[1] = 0.7;
        b[2] = 0.9;
        b[3] = 6.2;
    }
    return b;
}


struct matrix get_matrix_obj(int places, int variant) {
    struct matrix m = {4, get_matrix_with_sqrts(places, variant), get_b_for_sqrts_case(variant)};
    return m;
}

void run_tests(struct method* method) {
    const int places[] = {
        2,
        4,
        6,
        10
    };

    const int variants[] = { VARIANT_3, VARIANT_16 };

    for (int i = 0; i < sizeof(variants) / sizeof(int); i++) {
        for (int j = 0; j < sizeof(places) / sizeof(int); j++) {
        char buf[256];
        snprintf(buf, sizeof(buf), "results/%s_sqrt%d_var%d.txt", method->name, places[j], variants[i]);

        FILE *f = fopen(buf, "w");

        struct matrix m = get_matrix_obj(places[j], variants[i]);
        fill_aux_fields(&m, get_matrix_with_sqrts(-1, variants[i]));

        fprintf(f, "Initial a:\n");
        print_2d_array_to_file(f, m.n, m.a);

        fprintf(f, "Initial b:\n");
        print_array_to_file(f, m.n, m.b);

        method->func(&m);

        fprintf(f, "\nResult x:\n");
        print_array_to_file(f, m.n, m.x);

        fprintf(f, "\nResidual c:\n");
        print_array_to_file(f,  m.n, m.c);

        fprintf(f, "Part a after calculation:\n");
        print_2d_array_to_file(f, m.n, m.a);

        fprintf(f, "Part b after calculation:\n");
        print_array_to_file(f,  m.n, m.b);


        free_matrix(&m);
        fclose(f);
    }
    }
}

void _round_test() {
    FILE *f = fopen("results/_round_test.txt", "w");

    double arr[] = {
        sqrt(2),
        _round(sqrt(2), 1),
        _round(sqrt(2), 2),
        _round(sqrt(2), 3),
        _round(sqrt(2), 4),
        _round(sqrt(2), 5),
        _round(sqrt(2), 6),
        _round(sqrt(2), 7),
        _round(sqrt(2), 8),
        _round(sqrt(2), 9),
        _round(sqrt(2), 10)
    };

    for (int i = 0; i < sizeof(arr) / sizeof(double); i++) {
        fprintf(f, "%.10f\n", arr[i]);
    }

    fclose(f);
}