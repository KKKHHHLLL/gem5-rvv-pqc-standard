#include <stdio.h>
#include "defs.h"
#include <stdint.h>
#include <stdbool.h>
typedef uint32_t T;
typedef uint64_t U;

void reverse(uint32_t* x)
{
	*x = (((*x & 0xaaaaaaaa) >> 1) | ((*x & 0x55555555) << 1));
	*x = (((*x & 0xcccccccc) >> 2) | ((*x & 0x33333333) << 2));
	*x = (((*x & 0xf0f0f0f0) >> 4) | ((*x & 0x0f0f0f0f) << 4));
	*x = (((*x & 0xff00ff00) >> 8) | ((*x & 0x00ff00ff) << 8));

	*x = ((*x >> 16) | (*x << 16));
}

void General_bitreverse_coeff(const T a_in[R_SIZE_PRIME], T res[R_SIZE_PRIME]) {
	//right shift a_in to the rightmost
	uint32_t shift_mount = SEW - (R_BITS & (SEW-1));
	for (int i = R_SIZE_PRIME - 1; i >= 0; i--) {
		T f1 = a_in[i];
		T f0 = i == 0 ? 0 : a_in[i - 1];
		res[i] = (((((U)f1) << SEW | (U)f0) << shift_mount) >> SEW) & CHUNK_MASK;
	}
	//do bitreverse on every element
	for (int i = 0; i < R_SIZE_PRIME; i++)
		reverse(&res[i]);
	//reverse element order
	int left = 0;
	int right = R_SIZE_PRIME - 1;
	T temp;
	while (left < right) {
		temp = res[right];
		res[right] = res[left];
		res[left] = temp;
		left++;
		right--;
	}
}


void General_compute_control_bits(int32_t* delta,const uint32_t step_size, const T f[R_SIZE_PRIME], const T g[R_SIZE_PRIME], bool c[2 * S]) {
	T f_ = f[0];
	T g_ = g[0];
	for (int i = 0; i < step_size; i++) {
		bool alpha = g_ & 1;
		bool swap = (*delta > 0 ? true : false) && alpha;
		c[2 * i] = swap;
		c[2 * i + 1] = alpha;
		*delta = swap ? -(*delta) + 1 : (*delta) + 1;
		T f_temp = f_;
		f_ = swap ? g_ : f_;
		g_ = alpha ? (g_ ^ f_temp) >> 1 : g_ >> 1;
	}
}


void General_update_fg_or_vw(const bool c[2 * S], const uint32_t step_size,  const T f1,  const T f0,  const T g1,  const T g0,  const bool is_updating_fg,  T* r0,  T* r1) {
	U f = ((U)f1) << SEW | (U)f0;
	U g = ((U)g1) << SEW | (U)g0;
	for (int i = 0; i < step_size; i++) {
		U f_temp = f;
		f = c[2 * i] ? g : f;
		g = c[2 * i + 1] ? g ^ f_temp : g;
		if (is_updating_fg)
			g >>= 1;
		else
			f <<= 1;
	}
	if (is_updating_fg) {
		*r0 = f & CHUNK_MASK;
		*r1 = g & CHUNK_MASK;
	}
	else {
		*r0 = (f >> SEW) & CHUNK_MASK;
		*r1 = (g >> SEW) & CHUNK_MASK;
	}
}


void General_poly_inverse( const uint32_t step_size,  const T g_in[R_SIZE_PRIME],  T g_out[R_SIZE_PRIME]) {
	T f[R_SIZE_PRIME] = { 0 };
	T g[R_SIZE_PRIME] = { 0 };
	T* v = g_out;//g_out should be initialized to 0 before passed to this function
	T w[R_SIZE_PRIME] = { 0 };
	w[0] = 1, f[0] = 1, f[R_SIZE_PRIME - 1] = 1<<(R_BITS&(SEW-1));
	//reverse g_in to get g
	General_bitreverse_coeff(g_in, g);
	int32_t delta = 1;
	uint32_t Tau = 2 * R_BITS - 1;
	bool c[2 * S] = { 0 };

	//define a temp array for test
	uint32_t temp[R_SIZE_PRIME] = { 0 };

	while (Tau >= step_size) {
		General_compute_control_bits(&delta, step_size, f, g, c);
		//do stripmining on the chunks
		//update f g
		uint32_t chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T f0 = f[chunk_solved + j];
				T f1 = chunk_solved + j + 1 >= R_SIZE_PRIME ? 0 : f[chunk_solved + j + 1];
				T g0 = g[chunk_solved + j];
				T g1 = chunk_solved + j + 1 >= R_SIZE_PRIME ? 0 : g[chunk_solved + j + 1];
				T r0, r1;
				General_update_fg_or_vw(c, step_size, f1, f0, g1, g0, true, &r0, &r1);
				f[chunk_solved + j] = r0;
				g[chunk_solved + j] = r1;
			}
			chunk_solved += len;
		}

		//temp code for inspecting f
		General_bitreverse_coeff(g,temp);

		//update v w
		chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T v0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : v[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T v1 = v[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T w0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : w[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T w1 = w[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T r0, r1;
				General_update_fg_or_vw(c, step_size, v1, v0, w1, w0, false, &r0, &r1);
				v[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r0;
				w[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r1;
			}
			chunk_solved += len;
		}
		Tau -= step_size;
	}
	if (Tau > 0) {
		General_compute_control_bits(&delta, Tau, f, g, c);
		//update v w
		uint32_t chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T v0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : v[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T v1 = v[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T w0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : w[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T w1 = w[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T r0, r1;
				General_update_fg_or_vw(c, Tau, v1, v0, w1, w0, false, &r0, &r1);
				v[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r0;
				w[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r1;
			}
			chunk_solved += len;
		}
	}
	//shift v one bit right
	T temp_v[R_SIZE_PRIME] = { 0 };
	for (int i = 0; i < R_SIZE_PRIME; i++) {
		T v0 = v[i];
		T v1 = i + 1 >= R_SIZE_PRIME ? 0 : v[i + 1];
		temp_v[i] = ((((U)v1) << SEW | (U)v0) >> 1) & CHUNK_MASK;
	}
	//reverse v
	General_bitreverse_coeff(temp_v, v);
}

int main()
{
    uint32_t hin[R_SIZE_PRIME];
    uint32_t hout[R_SIZE_PRIME];
    printf("R=%d\n",R_BITS);
	unsigned long start, end;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
	General_poly_inverse(32,hin,hout);
	__asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("time spent:%lu\n",end-start);
    return 0;
}