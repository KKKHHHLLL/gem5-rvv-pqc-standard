#include "types.h"
#include <stdio.h>
#include <assert.h>
void gf2x_mul_base_port(OUT uint64_t *c,
                        IN const uint64_t *a,
                        IN const uint64_t *b)
{
  uint64_t       h = 0, l = 0, g1, g2, u[8];
  const uint64_t w  = 64;
  const uint64_t s  = 3;
  const uint64_t a0 = a[0];
  const uint64_t b0 = b[0];

  // Multiplying 64 bits by 7 can results in an overflow of 3 bits.
  // Therefore, these bits are masked out, and are treated in step 3.
  const uint64_t b0m = b0 & MASK(61);

  // Step 1: Calculate a multiplication table with 8 entries.
  u[0] = 0;
  u[1] = b0m;
  u[2] = u[1] << 1;
  u[3] = u[2] ^ b0m;
  u[4] = u[2] << 1;
  u[5] = u[4] ^ b0m;
  u[6] = u[3] << 1;
  u[7] = u[6] ^ b0m;

  // Step 2: Multiply two elements in parallel in positions i, i+s
  l = u[LSB3(a0)] ^ (u[LSB3(a0 >> 3)] << 3);
  h = (u[LSB3(a0 >> 3)] >> 61);

  for(size_t i = (2 * s); i < w; i += (2 * s)) {
    const size_t i2 = (i + s);

    g1 = u[LSB3(a0 >> i)];
    g2 = u[LSB3(a0 >> i2)];

    l ^= (g1 << i) ^ (g2 << i2);
    h ^= (g1 >> (w - i)) ^ (g2 >> (w - i2));
  }

  // Step 3: Multiply the last three bits.
  for(size_t i = 61; i < 64; i++) {
    uint64_t mask = (-((b0 >> i) & 1));
    l ^= ((a0 << i) & mask);
    h ^= ((a0 >> (w - i)) & mask);
  }

  c[0] = l;
  c[1] = h;
}

void karatzuba_add1_port(OUT uint64_t *alah,
                         OUT uint64_t *blbh,
                         IN const uint64_t *a,
                         IN const uint64_t *b,
                         IN const size_t    qwords_len)
{
  assert(qwords_len % REG_QWORDS == 0);

  REG_T va0, va1, vb0, vb1;

  for(size_t i = 0; i < qwords_len; i += REG_QWORDS) {
    va0 = LOAD(&a[i]);
    va1 = LOAD(&a[i + qwords_len]);
    vb0 = LOAD(&b[i]);
    vb1 = LOAD(&b[i + qwords_len]);

    STORE(&alah[i], va0 ^ va1);
    STORE(&blbh[i], vb0 ^ vb1);
  }
}

void karatzuba_add2_port(OUT uint64_t *z,
                         IN const uint64_t *x,
                         IN const uint64_t *y,
                         IN const size_t    qwords_len)
{
  assert(qwords_len % REG_QWORDS == 0);

  REG_T vx, vy;

  for(size_t i = 0; i < qwords_len; i += REG_QWORDS) {
    vx = LOAD(&x[i]);
    vy = LOAD(&y[i]);

    STORE(&z[i], vx ^ vy);
  }
}

void karatzuba_add3_port(OUT uint64_t *c,
                         IN const uint64_t *mid,
                         IN const size_t    qwords_len)
{
  assert(qwords_len % REG_QWORDS == 0);

  REG_T vr0, vr1, vr2, vr3, vt;

  uint64_t *c0 = c;
  uint64_t *c1 = &c[qwords_len];
  uint64_t *c2 = &c[2 * qwords_len];
  uint64_t *c3 = &c[3 * qwords_len];

  for(size_t i = 0; i < qwords_len; i += REG_QWORDS) {
    vr0 = LOAD(&c0[i]);
    vr1 = LOAD(&c1[i]);
    vr2 = LOAD(&c2[i]);
    vr3 = LOAD(&c3[i]);
    vt  = LOAD(&mid[i]);

    STORE(&c1[i], vt ^ vr0 ^ vr1);
    STORE(&c2[i], vt ^ vr2 ^ vr3);
  }
}

void karatzuba(OUT uint64_t *c,
                        IN const uint64_t *a,
                        IN const uint64_t *b,
                        IN const size_t    qwords_len,
                        IN const size_t    qwords_len_pad,
                        uint64_t *         sec_buf)
{
	//printf("len:%lu\n",qwords_len);
  if(qwords_len <= 1) {
    gf2x_mul_base_port(c, a, b);
    return;
  }

  const size_t half_qw_len = qwords_len_pad >> 1;

  // Split a and b into low and high parts of size n_padded/2
  const uint64_t *a_lo = a;
  const uint64_t *b_lo = b;
  const uint64_t *a_hi = &a[half_qw_len];
  const uint64_t *b_hi = &b[half_qw_len];

  // Split c into 4 parts of size n_padded/2 (the last ptr is not needed)
  uint64_t *c0 = c;
  uint64_t *c1 = &c[half_qw_len];
  uint64_t *c2 = &c[half_qw_len * 2];

  // Allocate 3 ptrs of size n_padded/2  on sec_buf
  uint64_t *alah = sec_buf;
  uint64_t *blbh = &sec_buf[half_qw_len];
  uint64_t *tmp  = &sec_buf[half_qw_len * 2];

  // Move sec_buf ptr to the first free location for the next recursion call
  sec_buf = &sec_buf[half_qw_len * 3];

  // Compute a_lo*b_lo and store the result in (c1|c0)
  karatzuba(c0, a_lo, b_lo, half_qw_len, half_qw_len, sec_buf);

  // If the real number of digits n is less or equal to n_padded/2 then:
  //     a_hi = 0 and b_hi = 0
  // and
  //     (a_hi|a_lo)*(b_hi|b_lo) = a_lo*b_lo
  // so we can skip the remaining two multiplications
  if(qwords_len > half_qw_len) {
    // Compute a_hi*b_hi and store the result in (c3|c2)
    karatzuba(c2, a_hi, b_hi, qwords_len - half_qw_len, half_qw_len, sec_buf);

    // Compute alah = (a_lo + a_hi) and blbh = (b_lo + b_hi)
    karatzuba_add1_port(alah, blbh, a, b, half_qw_len);

    // Compute (c1 + c2) and store the result in tmp
    karatzuba_add2_port(tmp, c1, c2, half_qw_len);

    // Compute alah*blbh and store the result in (c2|c1)
    karatzuba(c1, alah, blbh, half_qw_len, half_qw_len, sec_buf);

    // Add (tmp|tmp) and (c3|c0) to (c2|c1)
    karatzuba_add3_port(c0, tmp, half_qw_len);
  }
}

void gf2x_red_port(OUT pad_r_t *c, IN const dbl_pad_r_t *a)
{
  const uint64_t *a64 = (const uint64_t *)a;
  uint64_t *      c64 = (uint64_t *)c;

  for(size_t i = 0; i < R_QWORDS; i += REG_QWORDS) {
    REG_T vt0 = LOAD(&a64[i]);
    REG_T vt1 = LOAD(&a64[i + R_QWORDS]);
    REG_T vt2 = LOAD(&a64[i + R_QWORDS - 1]);

    vt1 = SLLI_I64(vt1, LAST_R_QWORD_TRAIL);
    vt2 = SRLI_I64(vt2, LAST_R_QWORD_LEAD);

    vt0 ^= (vt1 | vt2);

    STORE(&c64[i], vt0);
  }

  c64[R_QWORDS - 1] &= LAST_R_QWORD_MASK;
}

void gf2x_mod_mul(OUT pad_r_t *c,
                           IN const pad_r_t *a,
                           IN const pad_r_t *b)
{
  dbl_pad_r_t t = {0};
  uint64_t secure_buffer[SECURE_BUFFER_QWORDS]={0};

  karatzuba((uint64_t *)&t, (const uint64_t *)a, (const uint64_t *)b, R_QWORDS,
            R_PADDED_QWORDS, secure_buffer);

  gf2x_red_port(c, &t);
  //printf("1\n");//this line spent about 15661 cycles in plct-gem5
}

int main(){
	printf("R=%d\n",R_BITS);
	pad_r_t t;
	pad_r_t g;
    unsigned long start, end;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
	gf2x_mod_mul(&t, &g, &t);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
	//cycle cnt
	printf("take %lu cycles\n",end-start);

	return 0;
}