#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "instrumentation.h"
#define PERFORMANCE

#define print(...) printf(__VA_ARGS__)

#define M 64
#define N 64

// TODO:
// Priority queue (binary heap)
// Free list
// Cashe alloc
// int16_t for int
// Float for double, calculate z as c * x and keep track of c

// degub ldb

/* TYPEDEFS */
typedef struct simplex_t simplex_t;
typedef struct node_t node_t;
typedef struct list_node_t list_node_t;

simplex_t simplex_list[1000];
node_t node_list[1000];

// FUNCTION HEADERS
double xsimplex(int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y, int16_t* var, int16_t h);
double simplex(int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y);

// STRUCTS
struct simplex_t {
    int16_t m;          /* Constraints. */
    int16_t n;          /* Decision variables. */
    int16_t* var;       /* var[n+m+1] 0..n-1 are nonbasic. */
    double** a;     /* A. */
    double* b;      /* b. */
    double* x;      /* x[n+1] x. */
    double* c;      /* c. */
    double y;       /* y. */
};

struct node_t {
    int16_t m;          /* Constraints. */
    int16_t n;          /* Decision variables. */
    int16_t k;          /* Parent branches on x_k. */
    int16_t h;          /* Branch on x_h. */
    double xh;      /* x_h. */
    double ak;      /* Parent a_k. */
    double bk;      /* Parent b_k. */
    double* min;    /* Lower bounds. */
    double* max;    /* Upper bounds. */
    double** a;     /* A. */
    double* b;      /* b. */
    double* x;      /* x. */
    double* c;      /* c. */
    double z;       /* z. */
};

struct list_node_t {
    list_node_t* next;
    list_node_t* prev;
    node_t* data;
};

//FREE FUNCTIONS
void free_simplex(simplex_t* simplex) {
    for (int16_t i = 0; i < simplex->m; i++) {
        free(simplex->a[i]);
    }
    free(simplex->a);
    free(simplex->b);
    free(simplex->c);
    free(simplex->x);
    free(simplex->var);
    free(simplex);
}

void free_node(node_t* node) {
    for (int16_t i = 0; i < node->m+1; i++) {
        free(node->a[i]);
    }
    free(node->a);
    free(node->b);
    free(node->c);
    free(node->x);
    free(node->min);
    free(node->max);
    free(node);
}

node_t* initial_node(int16_t m, int16_t n, double** a, double* b, double* c) {
    struct node_t* p = calloc(1, sizeof(struct node_t));
    p->a = calloc(m+1, sizeof(double*));
    for (int16_t i = 0; i < m+1; i++) {
        p->a[i] = calloc(n+1, sizeof(double));
    }

    p->b = calloc(m+1, sizeof(double));
    p->c = calloc(n+1, sizeof(double));
    p->x = calloc(n+1, sizeof(double));
    p->min = calloc(n, sizeof(double));
    p->max = calloc(n, sizeof(double));
    p->m = m;
    p->n = n;

    for (int16_t i = 0; i < m; i++) {
        for (int16_t j = 0; j < n; j++) {
            p->a[i][j] = a[i][j];
        }
        p->b[i] = b[i];
    }

    for (int16_t i = 0; i < n; i++) {
        p->c[i] = c[i];
        p->min[i] = -INFINITY;
        p->max[i] = INFINITY;
    }

    return p;
}

node_t* extend(node_t* p, int16_t m, int16_t n, double** a, double* b, double* c, int16_t k, double ak, double bk) {
    struct node_t* q = calloc(1, sizeof(struct node_t));
    int16_t i, j;
    q->k = k;
    q->ak = ak;
    q->bk = bk;

    if (ak > 0 && p->max[k] < INFINITY) {
        q->m = p->m;
    } else if (ak < 0 && p->min[k] > 0) {
        q->m = p->m;
    } else {
        q->m = p->m + 1;
    }

    q->n = p->n;
    q->h = -1;

    q->a = calloc(q->m+1, sizeof(double*));
    for (int16_t i = 0; i < q->m+1; i++) {
        q->a[i] = calloc(q->n+1, sizeof(double));
    }

    q->b = calloc(q->m+1, sizeof(double));
    q->c = calloc(q->n+1, sizeof(double));
    q->x = calloc(q->n+1, sizeof(double));
    q->min = calloc(n, sizeof(double));
    q->max = calloc(n, sizeof(double));

    for (i = 0; i < p->n; i++) {
        q->min[i] = p->min[i];
        q->max[i] = p->max[i];
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            q->a[i][j] = a[i][j];
        }
        q->b[i] = b[i];
    }

    for (i = 0; i < n; i++) {
        q->c[i] = c[i];
    }

    if (ak > 0) {
        if (q->max[k] == INFINITY || bk < q->max[k]) {
            q->max[k] = bk;
        }
    } else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
        q->min[k] = -bk;
    }

    for (i = m, j = 0; j < n; j++) {
        if (q->min[j] > -INFINITY) {
            q->a[i][j] = -1;
            q->b[i] = -q->min[j];
            i += 1;
        }
        if (q->max[j] < INFINITY) {
            q->a[i][j] = 1;
            q->b[i] = q->max[j];
            i += 1;
        }
    }
    
    return q;
}

int is_integer(double* xp) {
    double x = *xp;
    double r = lround(x);
    if (fabs(r - x) < 1e-6) {
        *xp = r;
        return 1;
    } else {
        return 0;
    }
}

int integer(node_t* p) {
    for (int16_t i = 0; i < p->n; i++) {
        if (!is_integer(&p->x[i])) {
            return 0;
        }
    }
    return 1;
}

void bound(node_t* p, list_node_t** h, double* zp, double* x) {
    if (p->z > *zp) {
        *zp = p->z;

        for (int16_t i = 0; i < p->n; i++) {
            p->x[i] = x[i];
        }

        list_node_t* h_q = (*h);
        while (h_q != NULL) {
            list_node_t* next = h_q->next;
            if (h_q->data->z < p->z) {
                if (h_q->prev != NULL) {
                    h_q->prev->next = h_q->next;
                } else {
                    (*h) = h_q->next;
                }

                if (h_q->next != NULL) {
                    h_q->next->prev = h_q->prev;
                }
                free(h_q->data->min);
                free(h_q->data->max);
                free(h_q->data);
                free(h_q);
            }
            h_q = next;
        }
    }
}

double double_fraction(double fraction) {
    double integer_part = 0;
    double fractional_part;
    integer_part = modf(fraction, &fractional_part);
    return integer_part;
}

int branch(node_t* q, double z) {
    double min, max;
    
    if (q->z < z) {
        return 0;
    }

    double integer_part = 0;
    int16_t best_h = 0;
    double best_fractional = INFINITY;
    for (int16_t h = 0; h < q->n; h++) {
        if (!is_integer(&q->x[h])) {
            // fabs(modf(q->x[h], &integer_part) - 0.5) < 0.1
            // printf("TEST: %lf %lf\n", q->x[h], modf(q->x[h], &integer_part));
            if (q->min[h] == -INFINITY) {
                min = 0;
            } else {
                min = q->min[h];
            }
            max = q->max[h];
            if (floor(q->x[h]) < min || ceil(q->x[h]) > max) {
                continue;
            }
            
            double temp = fabs(modf(q->x[h], &integer_part) - 0.63);
            if (best_fractional > temp) {
                best_h = h;
                best_fractional = temp;
            }
        }
    }
    if (best_fractional != INFINITY) {
        q->h = best_h;
        q->xh = q->x[best_h];

        for (int16_t i = 0; i < q->m+1; i++) {
            free(q->a[i]);
        }
        free(q->a);
        free(q->b);
        free(q->c);
        free(q->x);
        return 1;
    }
    
    return 0;
}

void succ(node_t* p, list_node_t** h, int16_t m, int16_t n, double** a, double* b, double* c, int16_t k, int16_t ak, int16_t bk, double* zp, double* x) {
    node_t* q = extend(p, m, n, a, b, c, k, ak, bk);
    if (q == NULL) {
        return;
    }

    q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);
    if (isfinite(q->z)) {
        if (integer(q)) {
            bound(q, h, zp, x);
        } else if (branch(q, *zp)) {
            list_node_t* temp_h = calloc(1, sizeof(struct list_node_t));
            temp_h->data = q;
            temp_h->next = *h;

            if (*h != NULL) {
                (*h)->prev = temp_h;
            }
            *h = temp_h;
            return;
        }
    }
    free_node(q);
}

int init(simplex_t* s, int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y, int16_t* var) {
    int16_t i, k;

    s->m = m;
    s->n = n;
    s->a = a;
    s->b = b;
    s->c = c;
    s->x = x;
    s->y = y;
    s->var = var;

    if(s->var == NULL){
        s->var = calloc(m+n+1, sizeof(int16_t));
        for (i=0; i < m + n; i++){
            s->var[i] = i;
        }
    }

    for (k = 0, i = 1; i < m; i++) {
        if (s->b[i] < s->b[k]) {
            k = i;
        }
    }

    return k;
}

int select_nonbasic(simplex_t s) {
    // for (int16_t i = 0; i < s.n; i++) {
    //     if (s.c[i] > 1e-6) {
    //         return i;
    //     }
    // }
    // return -1;

    // Max c[i]:
    int16_t idx_max = 0;
    for (int16_t i = 1; i < s.n; i++) {
        if (s.c[i] > 1e-6 && s.c[i] > s.c[idx_max]) {
            idx_max = i;
        }
    }
    if (idx_max != 0 || s.c[0] > 1e-6)
        return idx_max;
    return -1;

    // Steepest edge
    // double max_edge = -1;
    // int16_t idx_max = 0;
    // for (int16_t i = 0; i < s.n; i++) {
    //     double sum = 0;
    //     for (int16_t j = 0; j < s.m; j++) {
    //         sum += pow(s.a[j][i], 2);
    //     }
    //     sum = s.c[i] / sqrt(sum + 1);
    //     if (sum > max_edge) {
    //         max_edge = sum;
    //         idx_max = i;
    //     }
    // }

    // if (max_edge > 0) {
    //     return idx_max;
    // }
    
    // return -1;
}

// inline __attribute__((always_inline))
void pivot(simplex_t* s, int16_t row, int16_t col) {
    //BENCHMARK_BEGIN(pivot);
    double** a = s->a;
    double* b = s->b;
    double* c = s->c;
    int16_t m = s->m;
    int16_t n = s->n;
    int16_t i, j, t;

    double a_rc = 1/a[row][col];

    // const int16_t FOLDS = 4;
    
    t = s->var[col];
    s->var[col] = s->var[n+row];
    s->var[n+row] = t;
    s->y = s->y + c[col] * b[row] * a_rc;

    //for (i = 0; i < col; i++) {
    //    c[i] = c[i] - c[col] * a[row][i] * a_rc;
    //}
    //for (i = col + 1; i < n; i++) {
    //    c[i] = c[i] - c[col] * a[row][i] * a_rc;
    //}
    double c_temp = c[col];
    #pragma GCC unroll(4)
    for (i = 0; i < n; i++) {
        c[i] = c[i] - c_temp * a[row][i] * a_rc;
    }
    c[col] = c_temp;
    c[col] = -c[col] * a_rc;

    
    // for (i = 0; i < row; i++) {
    //     b[i] = b[i] - a[i][col] * b[row] * a_rc;
    // }
    // for (i = row + 1; i < m; i++) {
    //     b[i] = b[i] - a[i][col] * b[row] * a_rc;
    // }
    double b_temp = b[row];
    #pragma GCC unroll(4)
    for (i = 0; i < m; i++) {
        b[i] = b[i] - a[i][col] * b_temp * a_rc;
    }
    b[row] = b_temp;

    //for (i = 0; i < row; i++) {
    //    for (j = 0; j < col; j++) {
    //        a[i][j] = a[i][j] - a[i][col] * a[row][j] * a_rc;
    //    }
    //    for (j = col + 1; j < n; j++) {
    //        a[i][j] = a[i][j] - a[i][col] * a[row][j] * a_rc;
    //    }
    //}
    //for (i = row + 1; i < m; i++) {
    //    for (j = 0; j < col; j++) {
    //        a[i][j] = a[i][j] - a[i][col] * a[row][j] * a_rc;
    //    }
    //    for (j = col + 1; j < n; j++) {
    //        a[i][j] = a[i][j] - a[i][col] * a[row][j] * a_rc;
    //    }
    //}
    for (i = 0; i < row; i++) {
        double a_temp = a[i][col];
        #pragma GCC unroll(4)
        for (j = 0; j < n; j++) {
            a[i][j] = a[i][j] - a_temp * a[row][j] * a_rc;
        }
        a[i][col] = a_temp;
    }
    for (i = row + 1; i < m; i++) {
        double a_temp = a[i][col];
        #pragma GCC unroll(4)
        for (j = 0; j < n; j++) {
            a[i][j] = a[i][j] - a_temp * a[row][j] * a_rc;
        }
        a[i][col] = a_temp;
    }

    double temp_a = a[row][col];
    for (i = 0; i < m; i++) {
        a[i][col] = -a[i][col] * a_rc;
    }
    a[row][col] = temp_a;

    for (i = 0; i < n; i++) {
        a[row][i] = a[row][i] * a_rc;
    }
    a[row][col] = temp_a;
    b[row] = b[row] * a_rc;
    a[row][col] = 1 * a_rc;
    //BENCHMARK_END_INTERVALL(pivot, 100000);
}

void prepare(simplex_t* s, int16_t k) {
    int16_t m = s->m;
    int16_t n = s->n;
    int16_t i;

    for (i = m+n; i > n; i--) {
        s->var[i] = s->var[i-1];
    }
    s->var[n] = m+n;
    n += 1;

    for (i = 0; i < m; i++) {
        s->a[i][n-1] = -1;
    }

    s->x = calloc(m+n, sizeof(double));
    s->c = calloc(n, sizeof(double));
    
    s->c[n-1] = -1;
    s->n = n;
    pivot(s, k, n-1);
}

int initial(simplex_t* s, int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y, int16_t* var) {
    int16_t i, j, k;
    double w;
    k = init(s, m, n, a, b, c, x, y, var);

    if (b[k] >= 0) {
        return 1; // check feasible
    }
    prepare(s, k);
    n = s->n;
    s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

    for (i = 0; i < m+n; i++) {
        if (s->var[i] == m+n-1) {
            if (fabs(s->x[i]) > 1e-6) {
                free(s->x);
                free(s->c);
                return 0;
            }
            else {
                break;
            }
        }
    }

    if (i >= n) {
        for (j = 0, k = 0; k < n; k++) {
            if (fabs(s->a[i-n][k]) > fabs(s->a[i-n][j])) {
                j = k;
            }
        }
        pivot(s, i-n, j);
        i = j;
    }

    if (i < n-1) {
        k = s->var[i];
        s->var[i] = s->var[n-1];
        s->var[n-1] = k;
        
        for (k = 0; k < m; k++) {
            w = s->a[k][n-1];
            s->a[k][n-1] = s->a[k][i];
            s->a[k][i] = w;
        }
    }

    free(s->c);
    s->c = c;
    s->y = y;

    for (k = n-1; k < n+m-1; k++) {
        s->var[k] = s->var[k+1];
    }

    s->n = s->n - 1;
    n = s->n;
    double* t = calloc(n, sizeof(double));
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            if (k == s->var[j]) {
                t[j] = t[j] + s->c[k];
                goto next_k;
            }
        }

        for (j = 0; j < m; j++) {
            if (s->var[n+j] == k) {
                break;
            }
        }

        s->y = s->y + s->c[k] * s->b[j];
        for (i = 0; i < n; i++) {
            t[i] = t[i] - s->c[k] * s->a[j][i];
        }
        
        next_k:;
    }

    for (i = 0; i < n; i++) {
        s->c[i] = t[i];
    }
    free(t);
    free(s->x);

    return 1;
}

double xsimplex(int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y, int16_t* var, int16_t h) {
    simplex_t s;
    int16_t i, row, col;

    if (!initial(&s, m, n, a, b, c, x, y, var)) {
        free(s.var);
        return NAN;
    }

    while ((col = select_nonbasic(s)) >= 0) {
        row = -1;
        double a_temp = 0;
        for (i = 0; i < m; i++) {
            if (a[i][col] > 1e-6 && (row < 0 || b[i] / a[i][col] < b[row] * a_temp)) {
                row = i;
                a_temp = 1/a[row][col];
            }
        }

        if (row < 0) {
            free(s.var);
            return INFINITY;
        }
        pivot(&s, row, col);
    }

    if (h == 0) {
        for (i = 0; i < n; i++) {
            if (s.var[i] < n) {
                x[s.var[i]] = 0;
            }
        }
        for (i = 0; i < m; i++) {
            if (s.var[n+i] < n) {
                x[s.var[n+i]] = s.b[i];
            }
        }
        free(s.var);
    } else {
        for (i = 0; i < n; i++) {
            x[i] = 0;
        }
        for (i = n; i < n+m; i++) {
            x[i] = s.b[i-n];
        }
    }
    return s.y;
}

double simplex(int16_t m, int16_t n, double** a, double* b, double* c, double* x, double y) {
    return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

double intopt(int16_t m, int16_t n, double** a, double* b, double* c, double* x) {
    node_t* p = initial_node(m, n, a, b, c);
    list_node_t* h = calloc(1, sizeof(struct list_node_t));
    h->data = p;
    double z = -INFINITY;
    
    p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);
    if (integer(p) || !isfinite(p->z)) {
        z = p->z;
        if (integer(p)) {
            for (int16_t i = 0; i < p->n; i++) {
                x[i] = p->x[i];
            }
            free_node(p);
            free(h);
            return z;
        }
    }
    branch(p, z);
    while (h != NULL) {
        p = h->data;
        list_node_t* old_node = h;
        if (h->next != NULL) {
            h->next->prev = NULL;
        }
        h = h->next;
        free(old_node);

        succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
        succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
        free(p->min);
        free(p->max);
        free(p);
    }
    if (z == -INFINITY) {
        return NAN;
    } else {
        return z;   
    }
}

// #define INCLUDE_MAIN
#ifdef INCLUDE_MAIN
int main(int argc, char const *argv[])
{
    int16_t m;
    int16_t n;

    scanf("%d %d\n", &m, &n);
    // printf("m = %d, n = %d\n", m, n);

    double** a;
    a = calloc(m, sizeof(double*));
    for (int16_t i = 0; i < m; i++) {
        a[i] = calloc(n+1, sizeof(double));
    }

    double* c = calloc(n, sizeof(double));
    double* b = calloc(m, sizeof(double));
    
    for (int16_t ic = 0; ic < n; ic++) {
        scanf("%lf", &c[ic]);
    }
    scanf("\n");

    for (int16_t i = 0; i < m; i++) {
        for (int16_t j = 0; j < n; j++) {
            scanf("%lf", &a[i][j]);
        }
        scanf("\n");
    }

    for (int16_t ib = 0; ib < m; ib++) {
        scanf("%lf", &b[ib]);
    }
    scanf("\n");

    // for (int16_t i = 0; i < m; i++) {
    //     for (int16_t j = 0; j < n; j++) {
    //         printf("%10.3lf", a[i][j]);
    //     }
    //     printf("\n");
    // }

    double* x = calloc(n+1, sizeof(double));
    double y = 0;

    double ans = intopt(m, n, a, b, c, x);

    //print_double_array(b, m);

    print("%.3lf\n", ans);

    // dealloc
    for (int16_t i = 0; i < m; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    // free(c);
    // free(x);

    return 0;
}
#endif
