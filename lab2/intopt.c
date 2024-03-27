#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* TYPEDEFS */
typedef struct simplex_t simplex_t;

#define PRINT_MATRIX(MATRIX, ROWS, COLS) do { \
    for (int i = 0; i < ROWS; i++) { \
        for (int j = 0; j < COLS; j++) { \
            printf("%d ", MATRIX[i][j]); \
        } \
        printf("\n"); \
    } \
} while(0)

struct simplex_t {
    int m;
    int n;
    int* var;
    double** a;
    double* b;
    double* x;
    double* c;
    double y;
};

int init(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    int i, k;
    
    *simplex = (simplex_t){.m = m, .n = n, .var = var, .a = a, .b = b, .x = x, .c = c, .y = y};

    if ((*simplex).var == NULL) {
        (*simplex).var = calloc(m+n+1, sizeof(int));

        for (i = 0; i < n+m; i++) {
            (*simplex).var[i] = i;
        }
    }

    for (k = 0, i = 1; i < m; i++) {
        if (b[i] < b[k]) {
            k = i;
        }
    }

    return k;
}

int select_nonbasic (simplex_t simplex) {
    for (int i = 0; i < simplex.n; i++) {
        if (simplex.c[i] > 1e-6) {
            return i;
        }
    }
    return -1;
}

int initial(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    init(simplex, m, n, a, b, c, x, y, var);
    return 1;
}

void pivot(simplex_t* simplex, int row, int col) {
    double** a = simplex->a;
    double* b = simplex->b;
    double* c = simplex->c;
    int m = simplex->m;
    int n = simplex->n;
    int i, j, t;

    t = (*simplex).var[col];
    (*simplex).var[col] = (*simplex).var[n+row];
    (*simplex).var[n+row] = t;
    (*simplex).y = (*simplex).y + c[col] * b[row] / a[row][col];

    for (i = 0; i < n; i++) {
        if(i != col) {
            c[i] = c[i] - c[col] * a[row][i] / a[row][col];
        }
    }
    c[col] = -c[col] / a[row][col];

    for (i = 0; i < m; i++) {
        if(i != row) {
            b[i] = b[i] - a[i][col] * b[row] / a[row][col];
        }
    }

    for(i = 0; i < m; i++) {
        if(i != row) {
            for(j = 0; j < n; j++) {
                if(j != col) {
                    a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
                }
            }
        }
    }

    for(i = 0; i < m; i++) {
        if (i != row) {
            a[i][col] = -a[i][col] / a[row][col];
        }
    }

    for(i = 0; i < n; i++) {
        if (i != col) {
            a[row][i] = a[row][i] / a[row][col];
        }
    }

    b[row] = b[row] / a[row][col];
    a[row][col] = 1 / a[row][col];
}

double glob;
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
    simplex_t s;
    int i, row, col;

    if (!initial(&s, m, n, a, b, c, x, y, var)) {
        free(s.var);
        return 0;
    }

    while ((col = select_nonbasic(s)) >= 0) {
        row = -1;
        for(i = 0; i < m; i++) {
            if(a[i][col] > 1e-6 && (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
                row = i;
            }
            glob++;
        }
        if (row < 0) {
            free(s.var);
            return INT64_MAX;
        }

        pivot(&s, row, col);
        glob = s.a[0][0];
    }

    if (h == 0) {
        for(i = 0; i < n; i++) {
            if (s.var[i] < n) {
                x[s.var[i]] = 0;
            }
        }
        for(i = 0; i < m; i++) {
            if (s.var[n+i] < n) {
                x[s.var[n+i]] = s.b[i];
            }
        }
        free(s.var);
    } else {
        for (int i = 0; i < n; i++)
            x[i] = 0;
        for (int i = n; i < n+m; i++)
            x[i] = s.b[i-n];
    }
    return s.y;
}

double simplex(int m, int n, double** a, double* b, double* c, double* x, double y) {
    return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

int main(int argc, char** argv)
{
    int m;
    int n;

    // Read size of matrix
    scanf("%d %d/n", &m, &n);
    //printf("m = %d, n = %d\n", m, n);

    // Allocate memory for pointers
    double** a;
    double* b;
    double* c;

    // Allocate memory
    a = calloc(m, sizeof(double*));
    for(int i = 0; i < m; i++)
        a[i] = calloc(n, sizeof(double));

    b = calloc(m, sizeof(double));
    c = calloc(n, sizeof(double));

    // Read from standard input
    for(int i = 0; i < n; i++)
        scanf("%lf", c+i);
    scanf("/n");

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++)
            scanf("%lf", a[i]+j);
        scanf("/n");
    }

    for(int i = 0; i < m; i++)
        scanf("%lf", b+i);
    scanf("/n");

    // Print values in "a"
    //for(int i = 0; i < m; i++) {
    //    for(int j = 0; j < n; j++)
    //        printf("%7.3lf ", a[i][j]);
    //    printf("\n");
    //}

    double* x = calloc(n+1, sizeof(double));
    double y;

    double awns = simplex(m, n, a, b, c, x, y);
    //printf("%lf\n", awns);

    // Free memory
    for(int i = 0; i < m; i++)
        free(a[i]);

    free(a);
    free(b);
    free(c);

    free(x);
    
    return 0;
}
