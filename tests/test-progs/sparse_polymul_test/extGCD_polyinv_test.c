#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <riscv_vector.h>
#include "defs.h"

void gen_hmatrix(int32_t* delta_ptr, uint64_t fm, uint64_t gm,uint32_t step_size,uint64_t h[4]) {
	for (int j = 1; j <= step_size; j++) {
		uint64_t g_msb = gm >> 63;
		if (g_msb == 0) {
			(*delta_ptr) += 1;
			gm <<= 1;
			h[2] <<= 1, h[3] <<= 1;
		}
		else {
			if ((*delta_ptr) > 0) {
				(*delta_ptr) = -(*delta_ptr) + 1;
				uint64_t temp = fm;
				fm = gm, gm = (temp ^ gm) << 1;
				uint64_t h0_temp = h[0];
				uint64_t h1_temp = h[1];
				h[0] = h[2], h[1] = h[3], h[2] = (h[2] ^ h0_temp) << 1, h[3] = (h[3] ^ h1_temp) << 1;
			}
			else {
				(*delta_ptr) += 1;
				gm = (gm ^ fm) << 1;
				h[2] = (h[2] ^ h[0]) << 1, h[3] = (h[3] ^ h[1]) << 1;
			}
		}
	}
}

//gin's first element is padded zero
void extGCD_polyinv(uint64_t gin[1+R_SIZE_PRIME_HALF],uint64_t gout[R_SIZE_PRIME_HALF]){
    uint8_t shiftcnt=63 - (R_BITS & 63);
    //initial some work variables
    uint64_t* f=(uint64_t*)calloc(R_SIZE_PRIME_HALF, sizeof(uint64_t));
    f[0] = 1 << shiftcnt, f[R_SIZE_PRIME_HALF - 1] = 1ULL << 63;
    uint64_t* g=(uint64_t*)calloc(R_SIZE_PRIME_HALF, sizeof(uint64_t));
    
    //perform csll on gin and put result into g
    uint64_t* vsrc1_addr=gin;
    uint64_t* vsrc2_addr=gin+2;
    uint64_t* vdst_addr=g+1;
    size_t vl=vsetvl_e64m1(2);//vl should be 2
    uint16_t temp;//has no meaningful usage
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(temp):[s]"r"(shiftcnt),[imm]"i"(0));
    for(int i=0;i+2<R_SIZE_PRIME_HALF;i+=2){
        vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,vl);
        vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,vl);
        vuint64m1_t vdst;
        __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(vsrc1));
        vse64_v_u64m1(vdst_addr,vdst,vl);
        //address update
        vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
    }
    g[0] = gin[1]<<shiftcnt;

    uint64_t* w=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    w[0]=1;
    uint64_t* v=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));

    int32_t delta = 0;
	uint32_t Tau = 2 * R_BITS - 1;

	uint32_t step_size = 64-1;

    //define some temporary stack spaces
    uint64_t* h0fv=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h1gw=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h2fv=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h3gw=(uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* buffer=(uint64_t*)calloc(R_SIZE_PRIME<<1, sizeof(uint64_t));

    while(Tau>=step_size){
        uint64_t h[4]={1,0,0,1};
        uint64_t gm=g[R_SIZE_PRIME_HALF-1];
        uint64_t fm=f[R_SIZE_PRIME_HALF-1];
        gen_hmatrix(&delta,fm,gm,step_size,h);

        //load h matrix into vector
        vuint64m1_t h0=vle64_v_u64m1(&h[0],1);
        vuint64m1_t h1=vle64_v_u64m1(&h[1],1);
        vuint64m1_t h2=vle64_v_u64m1(&h[2],1);
        vuint64m1_t h3=vle64_v_u64m1(&h[3],1);
        
        //do carry-less multiplication
        //1.h0 multiply f
        vsrc2_addr=f;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME_HALF;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h0));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<R_SIZE_PRIME;i+=2)h0fv[i>>1]=buffer[i]^buffer[i-1];
        // h0fv[0]=buffer[0];

        //2.h1 multiply g
        vsrc2_addr=g;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME_HALF;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h1));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<R_SIZE_PRIME;i+=2)h1gw[i>>1]=buffer[i]^buffer[i-1];
        // h1gw[0]=buffer[0];

        //3.h2 multiply f
        vsrc2_addr=f;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME_HALF;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h2));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<R_SIZE_PRIME;i+=2)h2fv[i>>1]=buffer[i]^buffer[i-1];
        // h2fv[0]=buffer[0];

        //4.h3 multiply g
        vsrc2_addr=g;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME_HALF;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h3));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<R_SIZE_PRIME;i+=2)h3gw[i>>1]=buffer[i]^buffer[i-1];
        // h3gw[0]=buffer[0];

        //update f, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=f;
        for(int i=0;i+2<R_SIZE_PRIME_HALF;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }
        f[R_SIZE_PRIME_HALF-1]=h0fv[R_SIZE_PRIME_HALF-1]^h1gw[R_SIZE_PRIME_HALF-1];

        //update g, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=g;
        for(int i=0;i+2<R_SIZE_PRIME_HALF;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }
        g[R_SIZE_PRIME_HALF-1]=h2fv[R_SIZE_PRIME_HALF-1]^h3gw[R_SIZE_PRIME_HALF-1];

        //do carry-less multiplication
        //1.h0 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h0));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h0fv[i>>1]=buffer[i]^buffer[i-1];
        // h0fv[0]=buffer[0];

        //2.h1 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h1));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h1gw[i>>1]=buffer[i]^buffer[i-1];
        // h1gw[0]=buffer[0];

        //3.h2 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h2));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h2fv[i>>1]=buffer[i]^buffer[i-1];
        // h2fv[0]=buffer[0];

        //4.h3 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h3));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h3gw[i>>1]=buffer[i]^buffer[i-1];
        // h3gw[0]=buffer[0];

        //update v, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=v;
        for(int i=0;i<R_SIZE_PRIME;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }
        
        //update w, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=w;
        for(int i=0;i<R_SIZE_PRIME;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }

        Tau-=step_size;
    }
    if(Tau>0){
        uint64_t h[4]={1,0,0,1};
        uint64_t gm=g[R_SIZE_PRIME_HALF-1];
        uint64_t fm=f[R_SIZE_PRIME_HALF-1];
        gen_hmatrix(&delta,fm,gm,Tau,h);

        //load h matrix into vector
        vuint64m1_t h0=vle64_v_u64m1(&h[0],1);
        vuint64m1_t h1=vle64_v_u64m1(&h[1],1);
        vuint64m1_t h2=vle64_v_u64m1(&h[2],1);
        vuint64m1_t h3=vle64_v_u64m1(&h[3],1);

        //do carry-less multiplication
        //1.h0 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h0));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h0fv[i>>1]=buffer[i]^buffer[i-1];
        // h0fv[0]=buffer[0];

        //2.h1 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h1));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h1gw[i>>1]=buffer[i]^buffer[i-1];
        // h1gw[0]=buffer[0];

        //3.h2 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h2));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h2fv[i>>1]=buffer[i]^buffer[i-1];
        // h2fv[0]=buffer[0];

        //4.h3 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        for(int i=0;i<R_SIZE_PRIME;i++){
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,1);
            vuint64m1_t vdst;
            __asm__ __volatile__ ("clmul64 %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(h3));
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc2_addr+=1,vdst_addr+=2;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        // for(int i=2;i<(R_SIZE_PRIME<<1);i+=2)h3gw[i>>1]=buffer[i]^buffer[i-1];
        // h3gw[0]=buffer[0];

        //update v, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=v;
        for(int i=0;i<R_SIZE_PRIME;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }
        
        //update w, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=w;
        for(int i=0;i<R_SIZE_PRIME;i+=2){
            vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
            vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
            vuint64m1_t vdst=vxor_vv_u64m1(vsrc1,vsrc2,2);
            vse64_v_u64m1(vdst_addr,vdst,2);
            //address update
            vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
        }
    }

    //shift v R_BITS right, process in equivalent left shift
    shiftcnt+=1;
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(temp):[s]"r"(shiftcnt),[imm]"i"(0));

    vsrc1_addr=&v[R_SIZE_PRIME_HALF-2];
    vsrc2_addr=&v[R_SIZE_PRIME_HALF];
    vdst_addr=gout;
    for(int i=0;i+2<R_SIZE_PRIME_HALF;i+=2){
        vuint64m1_t vsrc1=vle64_v_u64m1(vsrc1_addr,2);
        vuint64m1_t vsrc2=vle64_v_u64m1(vsrc2_addr,2);
        vuint64m1_t vdst;
        __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(vsrc1));
        vse64_v_u64m1(vdst_addr,vdst,2);
        //address update
        vsrc1_addr+=2,vsrc2_addr+=2,vdst_addr+=2;
    }
    uint64_t temp_array[2]={v[R_SIZE_PRIME-1],0};
    vuint64m1_t vsrc1=vle64_v_u64m1(&v[R_SIZE_PRIME-3],2);
    vuint64m1_t vsrc2=vle64_v_u64m1(temp_array,2);
    vuint64m1_t vdst;
    __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(vdst):[vt]"vr"(vsrc2),[vs]"vr"(vsrc1));
    vse64_v_u64m1(temp_array,vdst,2);
    gout[R_SIZE_PRIME_HALF-1]=temp_array[0];

    //free the stack spaces
    free(f),free(g),free(v),free(w);
    free(h0fv),free(h1gw),free(h2fv),free(h3gw),free(buffer);
}

typedef struct poly_buffer_in{
    union{
        uint32_t poly_in32[R_SIZE_PRIME+2];
        uint64_t poly_in64[1+R_SIZE_PRIME_HALF];
    };
}poly_buffer_in_t;

typedef struct poly_buffer_out{
    union{
        uint32_t poly_in32[R_SIZE_PRIME];
        uint64_t poly_in64[R_SIZE_PRIME_HALF];
    };
}poly_buffer_out_t;

//Only Simulate Arm NEON work procedure
int main(){
    poly_buffer_in_t h0;;
    poly_buffer_out_t inv_h0={{0}};
    unsigned long start, end;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    extGCD_polyinv(h0.poly_in64,inv_h0.poly_in64);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("total time spent by poly_inv_RVV:%lu\n",end-start);
    //print out 
    //for(int i=0;i<R_SIZE_PRIME;i++)printf("inv_h0[%d]=%u\n",i,inv_h0.poly_in32[i]);
}