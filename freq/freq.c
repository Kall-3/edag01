#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct word_t {
    char* word;
    int freq;
} word_t;

static inline __attribute((always_inline)) int isPrime(int num) {
    int i;
    if (num == 1)
        return 0;

    for (i = 2; i <= num / 2; i++) {
        if (num % i == 0)
            return 0;
    }

    return 1;
}

static inline __attribute__((always_inline)) int freq() {
    char* line = calloc(1048, sizeof(char));

    unsigned int used = 0, i, found = 0, line_nbr = 1;

    word_t* words = calloc(10000, sizeof(word_t));
    for (i = 0; i < 10000; i++)
        words[i].word = calloc(1048, sizeof(char));

    while (fgets(line, 1048, stdin) != NULL) {
        line[strlen(line)-1] = '\0';

        for (i = 0; i < used; i++) {
            if (!strcmp(words[i].word, line)) {
                if (words[i].freq == 0) {
                    break;
                }
                found = 1;
                break;
            }
        }
        if (!found) {
            if (!isPrime(line_nbr)) {
                strcpy(words[i].word, line);
                words[i].freq = 1;
                used++;
                printf("added %s\n", line);
            } else
                printf("trying to delete %s: not found\n", line);

        } else if (isPrime(line_nbr)) {
            words[i].freq = 0;
            printf("trying to delete %s: deleted\n", line);
        } else {
            printf("counted %s\n", line);
            words[i].freq++;
        }
        found = 0;
        line_nbr++;
    }

    word_t best = words[0];
    for (i = 1; i < used; i++) {
        if (words[i].freq > best.freq)
            best = words[i];
    }
    printf("result: %s %d\n", best.word, best.freq);

    free(line);
    for (i = 0; i < 10000; i++)
        free(words[i].word);
    free(words);

    return 0;
}

int main() {
    return freq();
}
