/* bigint.h — arbitrary-precision signed integers, from scratch.
 *
 * sign-magnitude representation; magnitude is an array of base-2^32
 * limbs, little-endian (limb 0 = least significant). Multiplication is
 * schoolbook below KARATSUBA_LIMBS and Karatsuba above. Division is
 * Knuth algorithm D. Decimal conversion walks base 10^9.
 */
#ifndef BIGINT_H
#define BIGINT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint32_t *limbs; /* little-endian, base 2^32 */
    size_t len;      /* used limbs, 0 for zero */
    size_t cap;
    int neg;         /* 1 negative; zero is always non-negative */
} bigint_t;

#define KARATSUBA_LIMBS 32

/* lifecycle */
void bi_init(bigint_t *x);
void bi_free(bigint_t *x);
void bi_copy(bigint_t *dst, const bigint_t *src);
void bi_reserve(bigint_t *x, size_t cap);

/* construction */
void bi_set_u64(bigint_t *x, uint64_t v);
void bi_set_i64(bigint_t *x, int64_t v);
/* parse decimal text with optional leading '-'; returns 0 on success */
int bi_from_decimal(bigint_t *x, const char *text);
/* always succeeds; *out is malloc'd, caller frees */
char *bi_to_decimal(const bigint_t *x);

/* comparison of magnitudes and full values */
int bi_cmp_mag(const bigint_t *a, const bigint_t *b); /* -1, 0, 1 */
int bi_cmp(const bigint_t *a, const bigint_t *b);     /* handles signs */

/* arithmetic (result may alias an input unless noted) */
void bi_add(bigint_t *out, const bigint_t *a, const bigint_t *b);
void bi_sub(bigint_t *out, const bigint_t *a, const bigint_t *b);
void bi_mul(bigint_t *out, const bigint_t *a, const bigint_t *b);
/* floored division semantics; divisor must be non-zero */
void bi_divmod(bigint_t *q, bigint_t *r, const bigint_t *u, const bigint_t *v);
void bi_mod_pow(bigint_t *out, const bigint_t *base, const bigint_t *exp, const bigint_t *mod);

/* helpers used by the CLI demos */
void bi_fact(bigint_t *out, unsigned n);
void bi_fib(bigint_t *out, unsigned n);

/* testing hooks (exposed for the test suite) */
void bi_trim(bigint_t *x);
size_t bi_bitlen(const bigint_t *x);

#endif /* BIGINT_H */
