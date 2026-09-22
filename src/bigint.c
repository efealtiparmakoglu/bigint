#include "bigint.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BASE 4294967296ULL /* 2^32 */

/* --- lifecycle ------------------------------------------------------------ */

void bi_init(bigint_t *x) {
    x->limbs = NULL;
    x->len = 0;
    x->cap = 0;
    x->neg = 0;
}

void bi_free(bigint_t *x) {
    free(x->limbs);
    x->limbs = NULL;
    x->len = 0;
    x->cap = 0;
    x->neg = 0;
}

void bi_reserve(bigint_t *x, size_t cap) {
    if (x->cap >= cap) return;
    size_t new_cap = x->cap ? x->cap : 8;
    while (new_cap < cap) new_cap *= 2;
    x->limbs = realloc(x->limbs, new_cap * sizeof(uint32_t));
    assert(x->limbs);
    x->cap = new_cap;
}

void bi_trim(bigint_t *x) {
    while (x->len > 0 && x->limbs[x->len - 1] == 0) x->len--;
    if (x->len == 0) x->neg = 0;
}

void bi_copy(bigint_t *dst, const bigint_t *src) {
    if (dst == src) return;
    bi_reserve(dst, src->len);
    memmove(dst->limbs, src->limbs, src->len * sizeof(uint32_t));
    dst->len = src->len;
    dst->neg = src->neg;
}

/* --- construction ---------------------------------------------------------- */

void bi_set_u64(bigint_t *x, uint64_t v) {
    bi_reserve(x, 3);
    x->limbs[0] = (uint32_t)(v & 0xFFFFFFFF);
    x->len = 1;
    if (v >> 32) {
        x->limbs[1] = (uint32_t)(v >> 32);
        x->len = 2;
    }
    x->neg = 0;
}

void bi_set_i64(bigint_t *x, int64_t v) {
    if (v >= 0) {
        bi_set_u64(x, (uint64_t)v);
        return;
    }
    bi_set_u64(x, (uint64_t)(-(v + 1)) + 1); /* avoid UB on INT64_MIN */
    x->neg = 1;
}

/* --- comparisons ------------------------------------------------------------ */

int bi_cmp_mag(const bigint_t *a, const bigint_t *b) {
    if (a->len != b->len) return a->len < b->len ? -1 : 1;
    for (size_t i = a->len; i-- > 0;) {
        if (a->limbs[i] != b->limbs[i]) return a->limbs[i] < b->limbs[i] ? -1 : 1;
    }
    return 0;
}

int bi_cmp(const bigint_t *a, const bigint_t *b) {
    if (a->neg != b->neg) return a->neg ? -1 : 1;
    int mag = bi_cmp_mag(a, b);
    return a->neg ? -mag : mag;
}

/* --- magnitude helpers ------------------------------------------------------ */

static void mag_add(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    const bigint_t *big = a->len >= b->len ? a : b;
    const bigint_t *small = big == a ? b : a;
    bi_reserve(out, big->len + 1);
    uint64_t carry = 0;
    size_t i = 0;
    for (; i < small->len; i++) {
        uint64_t sum = (uint64_t)big->limbs[i] + small->limbs[i] + carry;
        out->limbs[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
    for (; i < big->len; i++) {
        uint64_t sum = (uint64_t)big->limbs[i] + carry;
        out->limbs[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
    if (carry) {
        out->limbs[i] = (uint32_t)carry;
        out->len = i + 1;
    } else {
        out->len = big->len;
    }
    out->neg = 0;
}

/* a must be >= b in magnitude */
static void mag_sub(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    bi_reserve(out, a->len);
    int64_t borrow = 0;
    for (size_t i = 0; i < a->len; i++) {
        int64_t diff = (int64_t)a->limbs[i] - (i < b->len ? b->limbs[i] : 0) - borrow;
        borrow = diff < 0;
        if (diff < 0) diff += BASE;
        out->limbs[i] = (uint32_t)diff;
    }
    out->len = a->len;
    bi_trim(out);
}

static void mag_mul_schoolbook(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    bi_reserve(out, a->len + b->len + 1);
    memset(out->limbs, 0, (a->len + b->len + 1) * sizeof(uint32_t));
    out->len = a->len + b->len;
    for (size_t i = 0; i < a->len; i++) {
        uint64_t carry = 0;
        uint64_t ai = a->limbs[i];
        for (size_t j = 0; j < b->len; j++) {
            uint64_t cur = (uint64_t)out->limbs[i + j] + ai * b->limbs[j] + carry;
            out->limbs[i + j] = (uint32_t)cur;
            carry = cur >> 32;
        }
        size_t k = i + b->len;
        while (carry) {
            uint64_t cur = (uint64_t)out->limbs[k] + carry;
            out->limbs[k] = (uint32_t)cur;
            carry = cur >> 32;
            k++;
        }
    }
    bi_trim(out);
}

static void mag_add_into(bigint_t *acc, size_t offset, uint32_t value) {
    uint64_t sum = (uint64_t)acc->limbs[offset] + value;
    acc->limbs[offset] = (uint32_t)sum;
    size_t i = offset + 1;
    uint64_t carry = sum >> 32;
    while (carry) {
        sum = (uint64_t)acc->limbs[i] + carry;
        acc->limbs[i] = (uint32_t)sum;
        carry = sum >> 32;
        i++;
    }
}

/* Karatsuba: out = |a| * |b|. Splits at `half` limbs; the three recursive
 * products are z0 = a0*b0, z1 = (a1+a0)(b1+b0), z2 = a1*b1, assembled
 * limb-wise so the carries between the shifted products add naturally. */
static void mag_mul_karatsuba(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    if (a->len < KARATSUBA_LIMBS || b->len < KARATSUBA_LIMBS) {
        mag_mul_schoolbook(out, a, b);
        return;
    }
    size_t max = a->len > b->len ? a->len : b->len;
    size_t half = (max + 1) / 2;

    bigint_t a1, a0, b1, b0;
    bi_init(&a1); bi_init(&a0); bi_init(&b1); bi_init(&b0);

    a0.len = a->len < half ? a->len : half;
    bi_reserve(&a0, a0.len);
    memcpy(a0.limbs, a->limbs, a0.len * sizeof(uint32_t));
    a1.len = a->len > half ? a->len - half : 0;
    if (a1.len) {
        bi_reserve(&a1, a1.len);
        memcpy(a1.limbs, a->limbs + half, a1.len * sizeof(uint32_t));
    }
    b0.len = b->len < half ? b->len : half;
    bi_reserve(&b0, b0.len);
    memcpy(b0.limbs, b->limbs, b0.len * sizeof(uint32_t));
    b1.len = b->len > half ? b->len - half : 0;
    if (b1.len) {
        bi_reserve(&b1, b1.len);
        memcpy(b1.limbs, b->limbs + half, b1.len * sizeof(uint32_t));
    }

    bigint_t z0, z1, z2, sum_a, sum_b;
    bi_init(&z0); bi_init(&z1); bi_init(&z2); bi_init(&sum_a); bi_init(&sum_b);
    mag_mul_karatsuba(&z0, &a0, &b0);
    mag_add(&sum_a, &a1, &a0);
    mag_add(&sum_b, &b1, &b0);
    mag_mul_karatsuba(&z1, &sum_a, &sum_b);
    mag_mul_karatsuba(&z2, &a1, &b1);

    /* z1 -= z0 + z2 (mathematically always non-negative) */
    bigint_t diff, tmp;
    bi_init(&diff);
    bi_init(&tmp);
    mag_add(&diff, &z0, &z2);
    mag_sub(&tmp, &z1, &diff);
    bi_copy(&z1, &tmp);

    bi_reserve(out, 4 * half + 2);
    memset(out->limbs, 0, (4 * half + 2) * sizeof(uint32_t));
    out->len = 4 * half + 1;
    for (size_t i = 0; i < z0.len; i++) mag_add_into(out, i, z0.limbs[i]);
    for (size_t i = 0; i < z1.len; i++) mag_add_into(out, i + half, z1.limbs[i]);
    for (size_t i = 0; i < z2.len; i++) mag_add_into(out, i + 2 * half, z2.limbs[i]);
    out->neg = 0;
    bi_trim(out);

    bi_free(&a1); bi_free(&a0); bi_free(&b1); bi_free(&b0);
    bi_free(&z0); bi_free(&z1); bi_free(&z2); bi_free(&sum_a); bi_free(&sum_b);
    bi_free(&diff); bi_free(&tmp);
}

void bi_mul(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    bigint_t mag;
    bi_init(&mag);
    mag_mul_karatsuba(&mag, a, b);
    bi_copy(out, &mag);
    out->neg = mag.len != 0 && (a->neg != b->neg);
    bi_trim(out);
    bi_free(&mag);
}

/* --- public add/sub with signs ---------------------------------------------- */

void bi_add(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    if (a->neg == b->neg) {
        mag_add(out, a, b);
        out->neg = a->neg;
        bi_trim(out);
        return;
    }
    int cmp = bi_cmp_mag(a, b);
    if (cmp == 0) {
        bi_set_u64(out, 0);
        return;
    }
    const bigint_t *big = cmp > 0 ? a : b;
    const bigint_t *small = cmp > 0 ? b : a;
    mag_sub(out, big, small);
    out->neg = big->neg;
    bi_trim(out);
}

void bi_sub(bigint_t *out, const bigint_t *a, const bigint_t *b) {
    bigint_t flipped;
    bi_init(&flipped);
    bi_copy(&flipped, b);
    flipped.neg = !flipped.neg;
    bi_add(out, a, &flipped);
    bi_free(&flipped);
}

/* --- division: Knuth algorithm D -------------------------------------------- */

/* x >>= bits (magnitude only) */
static void mag_shr(bigint_t *x, unsigned bits) {
    if (bits == 0 || x->len == 0) return;
    size_t words = bits / 32;
    unsigned shift = bits % 32;
    if (words >= x->len) {
        x->len = 0;
        return;
    }
    if (shift == 0) {
        memmove(x->limbs, x->limbs + words, (x->len - words) * sizeof(uint32_t));
        x->len -= words;
        return;
    }
    for (size_t i = 0; i + words < x->len; i++) {
        uint64_t value = (uint64_t)x->limbs[i + words] >> shift;
        if (i + words + 1 < x->len) {
            value |= (uint64_t)x->limbs[i + words + 1] << (32 - shift);
        }
        x->limbs[i] = (uint32_t)value;
    }
    x->len -= words;
    bi_trim(x);
}

static void mag_mul_limb(bigint_t *out, const bigint_t *a, uint32_t m) {
    if (m == 0) {
        bi_set_u64(out, 0);
        return;
    }
    bi_reserve(out, a->len + 1);
    uint64_t carry = 0;
    for (size_t i = 0; i < a->len; i++) {
        uint64_t cur = (uint64_t)a->limbs[i] * m + carry;
        out->limbs[i] = (uint32_t)cur;
        carry = cur >> 32;
    }
    out->len = a->len;
    if (carry) {
        out->limbs[out->len++] = (uint32_t)carry;
    }
}

void bi_divmod(bigint_t *q_out, bigint_t *r_out, const bigint_t *u, const bigint_t *v) {
    assert(v->len > 0 && "division by zero");
    bigint_t q, r;
    bi_init(&q);
    bi_init(&r);

    int cmp = bi_cmp_mag(u, v);
    if (cmp < 0) {
        bi_set_u64(&q, 0);
        bi_copy(&r, u);
        r.neg = u->neg;
        bi_trim(&r);
    } else if (cmp == 0) {
        bi_set_u64(&q, 1);
        bi_set_u64(&r, 0);
    } else if (v->len == 1) {
        /* single-limb divisor: simple long division */
        uint64_t divisor = v->limbs[0];
        bi_reserve(&q, u->len);
        q.len = u->len;
        uint64_t rem = 0;
        for (size_t i = u->len; i-- > 0;) {
            uint64_t cur = (rem << 32) | u->limbs[i];
            q.limbs[i] = (uint32_t)(cur / divisor);
            rem = cur % divisor;
        }
        bi_trim(&q);
        bi_set_u64(&r, rem);
    } else {
        /* normalize: shift so the divisor's top limb has its high bit set */
        unsigned shift = 0;
        uint32_t top = v->limbs[v->len - 1];
        while ((top & 0x80000000u) == 0) {
            top <<= 1;
            shift++;
        }
        bigint_t vn, un;
        bi_init(&vn);
        bi_init(&un);
        mag_mul_limb(&vn, v, (uint32_t)1 << shift);
        mag_mul_limb(&un, u, (uint32_t)1 << shift);
        /* Knuth D needs one extra high limb on the dividend */
        bi_reserve(&un, un.len + 1);
        un.limbs[un.len] = 0;
        un.len += 1;

        q.limbs = calloc(un.len - vn.len + 2, sizeof(uint32_t));
        q.cap = un.len - vn.len + 2;

        for (size_t j = un.len - vn.len; j-- > 0;) {
            uint64_t numerator =
                ((uint64_t)un.limbs[j + vn.len] << 32) | un.limbs[j + vn.len - 1];
            uint64_t qhat = numerator / vn.limbs[vn.len - 1];
            uint64_t rhat = numerator % vn.limbs[vn.len - 1];
            while (qhat >= BASE ||
                   qhat * vn.limbs[vn.len - 2] > ((rhat << 32) | un.limbs[j + vn.len - 2])) {
                qhat--;
                rhat += vn.limbs[vn.len - 1];
                if (rhat >= BASE) break;
            }
            /* un[j..] -= qhat * vn */
            uint64_t borrow = 0, carry = 0;
            for (size_t i = 0; i < vn.len; i++) {
                uint64_t p = qhat * vn.limbs[i] + carry;
                carry = p >> 32;
                uint64_t sub = (uint64_t)un.limbs[i + j] - (uint32_t)p - borrow;
                if ((int64_t)sub < 0) {
                    sub += BASE;
                    borrow = 1;
                } else {
                    borrow = 0;
                }
                un.limbs[i + j] = (uint32_t)sub;
            }
            uint64_t sub = (uint64_t)un.limbs[j + vn.len] - carry - borrow;
            int needed_addback = 0;
            if ((int64_t)sub < 0) {
                sub += BASE;
                needed_addback = 1;
            }
            un.limbs[j + vn.len] = (uint32_t)sub;

            if (needed_addback) {
                qhat--;
                uint64_t c = 0;
                for (size_t i = 0; i < vn.len; i++) {
                    uint64_t sum = (uint64_t)un.limbs[i + j] + vn.limbs[i] + c;
                    un.limbs[i + j] = (uint32_t)sum;
                    c = sum >> 32;
                }
                un.limbs[j + vn.len] = (uint32_t)((uint64_t)un.limbs[j + vn.len] + c);
            }
            q.limbs[j] = (uint32_t)qhat;
        }

        q.len = un.len - vn.len;
        bi_trim(&q);
        r.limbs = calloc(vn.len + 1, sizeof(uint32_t));
        r.cap = vn.len + 1;
        memcpy(r.limbs, un.limbs, vn.len * sizeof(uint32_t));
        r.len = vn.len;
        mag_shr(&r, shift);
        bi_free(&vn);
        bi_free(&un);
    }

    q.neg = q.len != 0 && (u->neg != v->neg);
    r.neg = r.len != 0 && u->neg;
    bi_copy(q_out, &q);
    bi_copy(r_out, &r);
    bi_free(&q);
    bi_free(&r);
}

/* --- modular exponentiation -------------------------------------------------- */

/* Square-and-multiply with reduction after every step. Inputs are never
 * aliased with the outputs: each step multiplies into scratch `t`, then
 * reduces t into the accumulator. */
void bi_mod_pow(bigint_t *out, const bigint_t *base, const bigint_t *exp, const bigint_t *mod) {
    bigint_t result, b, e, q, r, t;
    bi_init(&result);
    bi_init(&b);
    bi_init(&e);
    bi_init(&q);
    bi_init(&r);
    bi_init(&t);

    bi_set_u64(&result, 1);
    bi_divmod(&q, &r, base, mod); /* b = base % mod */
    bi_copy(&b, &r);
    bi_copy(&e, exp);

    while (e.len > 0) {
        if (e.limbs[0] & 1) {
            bi_mul(&t, &result, &b);
            bi_divmod(&q, &result, &t, mod);
        }
        /* e >>= 1 */
        for (size_t i = 0; i < e.len; i++) {
            uint64_t value = e.limbs[i] >> 1;
            if (i + 1 < e.len) value |= (uint64_t)(e.limbs[i + 1] & 1) << 31;
            e.limbs[i] = (uint32_t)value;
        }
        bi_trim(&e);
        if (e.len == 0) break;
        bi_mul(&t, &b, &b);
        bi_divmod(&q, &r, &t, mod);
        bi_copy(&b, &r);
    }

    bi_copy(out, &result);
    bi_free(&result);
    bi_free(&b);
    bi_free(&e);
    bi_free(&q);
    bi_free(&r);
    bi_free(&t);
}

/* --- decimal ------------------------------------------------------------------ */

int bi_from_decimal(bigint_t *x, const char *text) {
    bi_set_u64(x, 0);
    int neg = 0;
    if (*text == '-') {
        neg = 1;
        text++;
    } else if (*text == '+') {
        text++;
    }
    if (!*text) return -1;
    for (const char *p = text; *p; p++) {
        if (!isdigit((unsigned char)*p)) return -1;
    }
    bigint_t ten, shifted, digit;
    bi_init(&ten);
    bi_init(&shifted);
    bi_init(&digit);
    bi_set_u64(&ten, 10);
    for (const char *p = text; *p; p++) {
        bi_set_u64(&digit, (uint64_t)(*p - '0'));
        bi_mul(&shifted, x, &ten);
        bi_add(x, &shifted, &digit);
    }
    bi_free(&ten);
    bi_free(&shifted);
    bi_free(&digit);
    x->neg = neg;
    bi_trim(x);
    return 0;
}

char *bi_to_decimal(const bigint_t *x) {
    if (x->len == 0) {
        char *zero = malloc(2);
        strcpy(zero, "0");
        return zero;
    }
    bigint_t mag, q, base, r;
    bi_init(&mag);
    bi_init(&q);
    bi_init(&base);
    bi_init(&r);
    bi_reserve(&mag, x->len);
    memcpy(mag.limbs, x->limbs, x->len * sizeof(uint32_t));
    mag.len = x->len;
    bi_set_u64(&base, 1000000000U);

    /* repeatedly divide by 10^9; the remainders are the base-10^9 digits,
     * least significant first */
    uint32_t *parts = malloc((x->len + 1) * sizeof(uint32_t));
    size_t count = 0;
    while (mag.len > 0) {
        bi_divmod(&q, &r, &mag, &base);
        parts[count++] = r.len ? r.limbs[0] : 0;
        bi_copy(&mag, &q);
        bi_trim(&mag);
    }

    char *buffer = malloc((size_t)count * 9 + 4);
    size_t used = (size_t)sprintf(buffer, "%u", parts[count - 1]);
    for (size_t i = count - 1; i-- > 0;) {
        used += (size_t)sprintf(buffer + used, "%09u", parts[i]);
    }
    if (x->neg) {
        memmove(buffer + 1, buffer, used + 1);
        buffer[0] = '-';
    }
    free(parts);
    bi_free(&mag);
    bi_free(&q);
    bi_free(&base);
    bi_free(&r);
    return buffer;
}

/* --- helpers ------------------------------------------------------------------- */

void bi_fact(bigint_t *out, unsigned n) {
    bi_set_u64(out, 1);
    bigint_t mul;
    bi_init(&mul);
    for (unsigned i = 2; i <= n; i++) {
        bi_set_u64(&mul, i);
        bi_mul(out, out, &mul);
    }
    bi_free(&mul);
}

void bi_fib(bigint_t *out, unsigned n) {
    bigint_t a, b, t;
    bi_init(&a);
    bi_init(&b);
    bi_init(&t);
    bi_set_u64(&a, 0);
    bi_set_u64(&b, 1);
    for (unsigned i = 0; i < n; i++) {
        bi_add(&t, &a, &b);
        bi_copy(&a, &b);
        bi_copy(&b, &t);
    }
    bi_copy(out, &a);
    bi_free(&a);
    bi_free(&b);
    bi_free(&t);
}

size_t bi_bitlen(const bigint_t *x) {
    if (x->len == 0) return 0;
    size_t bits = (x->len - 1) * 32;
    uint32_t top = x->limbs[x->len - 1];
    while (top) {
        bits++;
        top >>= 1;
    }
    return bits;
}
