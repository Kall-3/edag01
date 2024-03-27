#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
    int m;
    int n;

    // Read size of matrix
    scanf("%d %d/n", &m, &n);
    printf("m = %d, n = %d\n", m, n);

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
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++)
            printf("%7.3lf ", a[i][j]);
        printf("\n");
    }

    // Free memory
    for(int i = 0; i < m; i++)
        free(a[i]);

    free(b);
    free(c);
    
    return 0;
}
