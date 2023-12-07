#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "error.h"
#include "poly.h"

//#define MOD 2048
//#define SIZE 2049
//#define MMOD 2047
#define MOD 10000004
#define SIZE 10000005
#define MMOD 10000003

typedef struct poly_t {
    int coef;
    int exp;
} poly_t;

inline __attribute__((always_inline)) void insertion_sort(poly_t* restrict p, int length) {
    int i, j;
    poly_t key;

    for (i = 1; i < length; i++) {
        key = p[i];
        j = i - 1;

        // Move elements of array[0..i-1], that are greater than key, to one position ahead
        // of their current position
        while (j >= 0 && p[j].exp < key.exp) {
            p[j + 1] = p[j];
            j = j - 1;
        }
        p[j + 1] = key;
    }
}

inline __attribute__((always_inline)) poly_t* new_poly_from_string(const char* pol) {
    poly_t* p = calloc(SIZE, sizeof(poly_t));

    const char *current = pol;
    char token[128]; // make smaller
    char next_char;
    int negative = 1;
    int i = 0;
    int coef, exp;

    while (sscanf(current, "%127[^+-]", token) == 1) {
        coef = 1;
        exp = 0;
        if (sscanf(token, "%dx^%d", &coef, &exp) == 2 || sscanf(token, "x^%d", &exp)) {
            // x^ term with no coef
        } else if (sscanf(token, "%d%c", &coef, &next_char) == 2 && next_char == 'x') {
            exp = 1;  // x term
        } else if (sscanf(token, "%d", &coef) == 1) {
            // constant term
        }

        // Add term to polynomial
        p[i].coef = coef * negative;
        p[i].exp = exp;
        negative = 1;

        // Move to the next term
        current += strlen(token);
        if (current[0] == '-') {
            negative = -1;
        }
        //current += (*current != '\0');
        current += !!(*current); // Skip the + or - sign
        // current++;  
        i++;
    }
    return p;
}

inline __attribute__((always_inline)) void free_poly(poly_t* restrict p) {
    free(p);
}

inline __attribute__((always_inline)) poly_t* mul(poly_t* restrict p, poly_t* restrict q) {
    poly_t* r = calloc(1024, sizeof(poly_t));
    int i, j, exp;

    for (i = 0; i <= 1024; i++) {
        for (j = 0; j <= 1024; j++) {
            exp     = p[i].exp + q[j].exp;
            
            r[exp & 1023].coef += p[i].coef * q[j].coef;
            r[exp & 1023].exp = exp;
        }
    }

    insertion_sort(r, 1024);
    return r;
}

inline __attribute__((always_inline)) void print_poly(poly_t* restrict p) {
    int a;

    for (int i = 0; p[i].coef != 0; i++) {
        if ((a = abs(p[i].coef)) > 1 || (a == 1 && p[i].exp == 0))
            printf("%d", a);

        if (p[i].exp > 1)
            printf("x^%d", p[i].exp);
        else if (p[i].exp == 1)
            printf("x");
        
        //dlm = p[i + 1].coef > 0 ? " + " : p[i + 1].coef < 0 ? " - " : "";
        if (p[i + 1].coef > 0) {
            printf(" + ");
        } else if (p[i + 1].coef < 0) {
            printf(" - ");
        }
    }
    printf("\n");
}
