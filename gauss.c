#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

struct matrix {
    int n;      // размерность матрицы
    double** a; // основная часть матрицы
    double* b;  // расширенная часть
    double* x;  // корни
    double* c;  // вектор невязки

    int* perm;        // вектор перестановок
    double** t;       // вспомогательная матрица
    double** a_exact; // значения a без округления
    double* b_copy;   // исходные значения b
};

void print_array(int n, double* a);
void print_2d_array(int n, double** a);
void print_array_to_file(FILE* f, int n, double* a);
void print_2d_array_to_file(FILE* f, int n, double** a);

struct matrix get_test_matrix();
void fill_aux_fields(struct matrix *m, double** matrix_without_rounding);
struct matrix get_matrix_from_input();
void free_matrix(struct matrix *m);

void calculate(struct matrix *m) {
    int n = m->n;
    double** a = m->a;
    double* b = m->b;
    double** t = m->t;
    double* x = m->x;

    for (int k = 0; k < n-1; k++) {
        for (int i = k+1; i < n; i++) {
            t[i][k] = a[i][k] / a[k][k];
            b[i] = b[i] - t[i][k] * b[k];
            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] - t[i][k] * a[k][j];
            }
        }
    }

    x[n-1] = b[n-1] / a[n-1][n-1];
    for (int k = n-2; k > -1; k--) {
        double s = 0.;
        for (int j = k+1; j < n; j++) {
            s = s + a[k][j] * x[j];
        }
        x[k] = (b[k] - s) / a[k][k];
    }

    double** a_exact = m->a_exact;
    double* b_copy = m->b_copy;
    double* c = m->c;

    for (int i = 0; i < n; i++) {
        double s = 0.;
        for (int j = 0; j < n; j++) {
            s += a_exact[i][j] * x[j];
        }

        c[i] = s - b_copy[i];
    }
}

void calculate_with_pivoting(struct matrix* m) {
    int n = m->n;
    double** a = m->a;
    double* b = m->b;
    double* x = m->x;

    for (int k = 0; k < n-1; k++) {
        double max = fabs(a[k][k]);
        int max_i = k;
        int max_j = k;
        for (int i = k; i < n; i++) {
            for (int j = k; j < n; j++) {
                if (fabs(a[i][j]) > max) {
                    max_i = i;
                    max_j = j;
                    max = fabs(a[i][j]);
                }
            }
        }

        double* swap_buf = malloc(sizeof(double) * n);
        swap_buf = a[max_i];
        a[max_i] = a[k];
        a[k] = swap_buf;

        double temp = b[k];
        b[k] = b[max_i];
        b[max_i] = temp;

        int temp_idx = m->perm[k];
        m->perm[k] = m->perm[max_j];
        m->perm[max_j] = temp_idx;

        for (int i = 0; i < n; i++) {
            double temp = a[i][k];
            a[i][k] = a[i][max_j];
            a[i][max_j] = temp;
        }

        for (int i = k+1; i < n; i++) {
            double t = a[i][k] / a[k][k];
            b[i] = b[i] - t * b[k];

            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] - t * a[k][j];
            }
        }

        swap_buf = NULL;
        free(swap_buf);
    }

    double* perm_buf = malloc(sizeof(double) * n);
    // if (fabs(a[n-1][n-1]) <= 1e-5 && fabs(b[n-1]) >= 1e-5) {
    //     printf("Решений нет!\n");
    //     exit(1);
    // }
    perm_buf[n-1] = b[n-1] / a[n-1][n-1];
    for (int k = n-2; k > -1; k--) {
        double s = 0.;
        for (int j = k+1; j < n; j++) {
            s = s + a[k][j] * perm_buf[j];
        }
        perm_buf[k] = (b[k] - s) / a[k][k];
    }

    for (int i = 0; i < n; i++) {
        x[m->perm[i]] = perm_buf[i];
    }
    perm_buf = NULL;
    free(perm_buf);

    double** a_exact = m->a_exact;
    double* b_copy = m->b_copy;
    double* c = m->c;

    for (int i = 0; i < n; i++) {
        double s = 0.;
        for (int j = 0; j < n; j++) {
            s += a_exact[i][j] * x[j];
        }

        c[i] = s - b_copy[i];
    }
}

void print_array_to_file(FILE* f, int n, double* a) {
    for (int i = 0; i < n; i++) {
        fprintf(f, "%14.10f ", a[i]);
    }
    fprintf(f, "\n");
}
void print_2d_array_to_file(FILE* f, int n, double** a) {
    for (int i = 0; i < n; i++) {
        print_array_to_file(f, n, a[i]);
    }
}

void print_array(int n, double* a) {
    for(int i = 0; i < n; i++) {
        printf("%f ", a[i]);
    }
    printf("\n");
}

void print_2d_array(int n, double** a) {
    for (int i = 0; i < n; i++) {
        print_array(n, a[i]);
    }
}

void fill_aux_fields(struct matrix *m, double** matrix_without_rounding) {
    size_t size_n = sizeof(double) * m->n;
    size_t size_ptr_n = sizeof(double*) * m->n;

    m->x = malloc(size_n);
    memset(m->x, 0, size_n);

    m->c = malloc(size_n);
    memset(m->c, 0, size_n);

    m->perm = malloc(size_n);
    for (int i = 0; i < m->n; i++) {
        m->perm[i] = i;
    }

    m->t = malloc(size_ptr_n);
    m->a_exact = malloc(size_ptr_n);
    for (int i = 0; i < m->n; i++) {
        m->t[i] = malloc(size_n);
        memset(m->t[i], 0, size_n);

        m->a_exact[i] = malloc(size_n);
        for (int j = 0; j < m->n; j++) {
            m->a_exact[i][j] = matrix_without_rounding[i][j];
        }
    }

    m->b_copy = malloc(size_n);
    memcpy(m->b_copy, m->b, size_n);
}

struct matrix get_matrix_from_input() {
    struct matrix m;

    scanf("%d", &m.n);

    size_t size_n = sizeof(double) * m.n;
    size_t size_ptr_n = sizeof(double*) * m.n;

    m.a = malloc(size_ptr_n);
    for (int i = 0; i < m.n; i++) {
        m.a[i] = malloc(size_n);

        for (int j = 0; j < m.n; j++) {
            scanf("%lf", &m.a[i][j]);
        }
    }
    m.b = malloc(size_n);
    for (int i = 0; i < m.n; i++) {
        scanf("%lf", &m.b[i]);
    }
    return m;
}

void free_matrix(struct matrix *m) {
    free(m->b);
    free(m->x);
    free(m->c);
    free(m->perm);

    for (int i = 0; i < m->n; i++) {
        free(m->a[i]);
        free(m->t[i]);
        free(m->a_exact[i]);
    }
    free(m->a);
    free(m->t);
    free(m->a_exact);

    free(m->b_copy);
}