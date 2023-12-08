#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include <assert.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : -(a))

#define BLK "\e[0;30m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"
#define YEL "\e[0;33m"
#define BLU "\e[0;34m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define WHT "\e[0;37m"

// Regular bold text
#define BBLK "\e[1;30m"
#define BRED "\e[1;31m"
#define BGRN "\e[1;32m"
#define BYEL "\e[1;33m"
#define BBLU "\e[1;34m"
#define BMAG "\e[1;35m"
#define BCYN "\e[1;36m"
#define BWHT "\e[1;37m"

// Regular underline text
#define UBLK "\e[4;30m"
#define URED "\e[4;31m"
#define UGRN "\e[4;32m"
#define UYEL "\e[4;33m"
#define UBLU "\e[4;34m"
#define UMAG "\e[4;35m"
#define UCYN "\e[4;36m"
#define UWHT "\e[4;37m"

// Regular background
#define BLKB "\e[40m"
#define REDB "\e[41m"
#define GRNB "\e[42m"
#define YELB "\e[43m"
#define BLUB "\e[44m"
#define MAGB "\e[45m"
#define CYNB "\e[46m"
#define WHTB "\e[47m"

// High intensty background
#define BLKHB "\e[0;100m"
#define REDHB "\e[0;101m"
#define GRNHB "\e[0;102m"
#define YELHB "\e[0;103m"
#define BLUHB "\e[0;104m"
#define MAGHB "\e[0;105m"
#define CYNHB "\e[0;106m"
#define WHTHB "\e[0;107m"

// High intensty text
#define HBLK "\e[0;90m"
#define HRED "\e[0;91m"
#define HGRN "\e[0;92m"
#define HYEL "\e[0;93m"
#define HBLU "\e[0;94m"
#define HMAG "\e[0;95m"
#define HCYN "\e[0;96m"
#define HWHT "\e[0;97m"

// Bold high intensity text
#define BHBLK "\e[1;90m"
#define BHRED "\e[1;91m"
#define BHGRN "\e[1;92m"
#define BHYEL "\e[1;93m"
#define BHBLU "\e[1;94m"
#define BHMAG "\e[1;95m"
#define BHCYN "\e[1;96m"
#define BHWHT "\e[1;97m"
#define RESET "\033[0m"

#define PERFORMANCE
#ifdef PERFORMANCE

#  define NS_PER_SECOND 1000000000

void get_delta_time(struct timespec t1, struct timespec t2, struct timespec* td)
{
  td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
  td->tv_sec = t2.tv_sec - t1.tv_sec;
  if (td->tv_sec > 0 && td->tv_nsec < 0) {
    td->tv_nsec += NS_PER_SECOND;
    td->tv_sec--;
  } else if (td->tv_sec < 0 && td->tv_nsec > 0) {
    td->tv_nsec -= NS_PER_SECOND;
    td->tv_sec++;
  }
}

#  define INTERVALL 1000

#  define BENCHMARK_BEGIN(name) \
    struct timespec benchmark_start_##name, benchmark_end_##name, \
        benchmark_delta_##name; \
    _Thread_local static int benchmark_count_##name = 0; \
    _Thread_local static double benchmark_total_##name = 0.0; \
    _Thread_local static double benchmark_min_##name = 10000.0; \
    _Thread_local static double benchmark_max_##name = 0.0; \
    clock_gettime(CLOCK_REALTIME, &benchmark_start_##name)

#  define PRINT_BENCHMARK(name) \
    printf("Benchmark: [" BGRN #name RESET "]: took " CYN "%.5f " RESET \
           "s " WHT \
           "(min: %.5f, max: %.5f, count: %d, total: %.5f, avg: %.5f)" RESET \
           "\n", \
           benchmark_tmp_##name, \
           benchmark_min_##name, \
           benchmark_max_##name, \
           benchmark_count_##name, \
           benchmark_total_##name, \
           benchmark_total_##name / benchmark_count_##name);

#  define BENCHMARK_END_HELPER(name) \
    clock_gettime(CLOCK_REALTIME, &benchmark_end_##name); \
    get_delta_time(benchmark_start_##name, \
                   benchmark_end_##name, \
                   &benchmark_delta_##name); \
    double benchmark_tmp_##name = \
        benchmark_delta_##name.tv_sec + benchmark_delta_##name.tv_nsec * 1e-9; \
    benchmark_count_##name++; \
    benchmark_total_##name += benchmark_tmp_##name; \
    benchmark_min_##name = MIN(benchmark_min_##name, benchmark_tmp_##name); \
    benchmark_max_##name = MAX(benchmark_max_##name, benchmark_tmp_##name);

#  define BENCHMARK_END(name) \
    BENCHMARK_END_HELPER(name); \
    PRINT_BENCHMARK(name)

#  define BENCHMARK_END_INTERVALL(name, intervall) \
    BENCHMARK_END_HELPER(name); \
    if (benchmark_count_##name % intervall == 0) { \
      PRINT_BENCHMARK(name); \
    }

#else
#  define BENCHMARK_BEGIN(...)
#  define BENCHMARK_END(...)
#  define BENCHMARK_END_INTERVALL(...)
#  define BENCHMARK_END_HELPER(...)
#  define PRINT_BENCHMARK(...)

#endif  // PERFORMANCE

#endif  // PERFORMANCE_H

