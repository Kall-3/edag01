#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// --------------------//
double sumArrayElements(double* array, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += array[i];
    }
    return sum;
}
//--------------------------//

/* TYPEDEFS */
typedef struct simplex_t simplex_t;

void printVector(double* vector, int length) {
    for (int i = 0; i < length; i++) {
        printf("%.2f ", vector[i]);
    }
    printf("\n\n");
}

void printMatrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

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

// Initialize functions
int init(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
int select_nonbasic (simplex_t simplex);
void prepare (simplex_t* simplex, int k);
int initial(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var);
void pivot(simplex_t* simplex, int row, int col);
double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h);
double simplex(int m, int n, double** a, double* b, double* c, double* x, double y);


int init(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    int i, k;

    *simplex = (simplex_t){.m = m, .n = n, .var = var, .a = a, .b = b, .x = x, .c = c, .y = y};

    if (simplex->var == NULL) {
        simplex->var = calloc(m+n+1, sizeof(int));

        for (i = 0; i < n+m; i++) {
            simplex->var[i] = i;
        }
    }

    for (k = 0, i = 1; i < m; i++) {
        if (simplex->b[i] < simplex->b[k]) {
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

void prepare(simplex_t* simplex, int k) {
    int m = simplex->m;
    int n = simplex->n;
    int i;

    for (i = m + n; i > n; i--) {
        simplex->var[i] = simplex->var[i-1];
    }
    simplex->var[n] = m + n;
    n = n + 1;

    printf("%d\n", n);
    // printMatrix(simplex->a, m, n-1);
    for (i = 0; i < m; i++) {
        simplex->a[i][n-1] = -1;
    }
    // printMatrix(simplex->a, m, n);

    free(simplex->x);
    free(simplex->c);

    simplex->x = calloc(m + n, sizeof(double));
    simplex->c = calloc(n, sizeof(double));

    simplex->c[n - 1] = -1;
    simplex->n = n;

    pivot(simplex, k, n-1);
    // printf("Prepare end\n");
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

int initial(simplex_t* simplex, int m, int n, double** a, double* b, double* c, double* x, double y, int* var) {
    int i, j, k;
    double w;

    k = init(simplex, m, n, a, b, c, x, y, var);

    // printVector(simplex->b, m);
    if (b[k] >= 0) {
        return 1;
    }
    prepare(simplex, k);
    n = simplex->n;
    simplex->y = xsimplex(m, n, simplex->a, simplex->b, simplex->c, simplex->x, 0, simplex->var, 1);

    for (i = 0; i < m+n; i++) {
        if (simplex->var[i] = m + n - 1) {
            if (abs(simplex->x[i]) > 1e-6) {
                free(simplex->x);
                free(simplex->c);
                return 0;
            } else {
                break;
            }
        }
    }

    if (i >= n) {
        for (j = 0, k = 0; k < n; k++) {
            if (abs(simplex->a[i-n][k]) > abs(simplex->a[i-n][j])) {
                j = k;
            }
        }
        pivot(simplex, i-n, j);
        i = j;
    }
    if (i < n-1) {
        k = simplex->var[i];
        simplex->var[i] = simplex->var[n-1];
        simplex->var[n-1] = k;
        for (k = 0; k < m; k++) {
            w = simplex->a[k][n-1];
            simplex->a[k][n-1] = simplex->a[k][i];
            simplex->a[k][i] = w;
        }
    }

    free(simplex->c);
    simplex->c = c;
    simplex->y = y;

    for (k = n-1; k < n+m-1; k++) {
        simplex->var[k] = simplex->var[k+1];
    }
    simplex->n = simplex->n-1;
    n = simplex->n;
    double* t = calloc(n, sizeof(double));

    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            if (k == simplex->var[j]) {
                t[j]+=simplex->c[k];
                goto next_k;
            }
        }
        for (j = 0; j < m; m++) {
            if (simplex->var[n+j] == k) {
                break;
            }
        }
        simplex->y = simplex->y + (simplex->c[k] * simplex->b[j]);

        for (i = 0; i < n; i++) {
            t[i] = t[i] - (simplex->c[k] * simplex->a[j][i]);
        }

    next_k:;
    }

    for (i = 0; i < n; i++) {
        simplex->c[i] = t[i];
    }
    free(t);
    free(simplex->x);

    return 1;
}

double xsimplex(int m, int n, double** a, double* b, double* c, double* x, double y, int* var, int h) {
    simplex_t s;
    int i, row, col;

    if (!initial(&s, m, n, a, b, c, x, y, var)) {
        free(s.var);
        return 0;
    }

    double* coords = calloc(m, sizeof(double)); // [0, 0]
    int* nonbasic_order = calloc(n, sizeof(int));
    int counter = 0;

    printf("First\n");
    printMatrix(a,m,n);
    printVector(b,m);
    printVector(c,n);
    printVector(coords, n);
    printf("%lf\n--------\n", s.y);

    while ((col = select_nonbasic(s)) >= 0 || sumArrayElements(coords, counter) < 5) {
        row = -1;
        for(i = 0; i < m; i++) {
            if(a[i][col] > 1e-6 && (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
                row = i;
            }
        }
        if (row < 0) {
            free(s.var);
            free(coords);
            free(nonbasic_order);
            return INT64_MAX;
        }

        pivot(&s, row, col);

        nonbasic_order[counter] = row;

        for (i = 0; i <= counter; i++) {
            coords[i] = s.b[nonbasic_order[i]]; // The nonbasic extraction gives one coordinate
        }

        counter++;

        printMatrix(a,m,n);
        printVector(b,m);
        printVector(c,n);
        printVector(coords, n);
        printf("%lf\n--------\n", s.y);
    }

    free(coords);
    free(nonbasic_order);

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
    printf("START\n");

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
    a = calloc(m+1, sizeof(double*));
    for(int i = 0; i < m; i++)
        a[i] = calloc(n+1, sizeof(double));

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
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n+1; j++)
            printf("%7.3lf ", a[i][j]);
        printf("\n");
    }
    // printMatrix(a, m, n);

    double* x = calloc(n+1, sizeof(double));
    double y = 0;

    double awns = simplex(m, n, a, b, c, x, y);
    // printf("%lf\n", awns);

    // Free memory
    for(int i = 0; i < m; i++)
        free(a[i]);

    free(a);
    free(b);
    free(c);

    free(x);
    
    return 0;
}
