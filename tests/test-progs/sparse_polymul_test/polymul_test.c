#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <riscv_vector.h>
#include "defs.h"

typedef struct poly_buffer{
    union{
        struct{
            uint32_t additional_word[ADDITIONAL];
            uint32_t val_word[R_SIZE_PRIME];
            uint32_t padded_word[PADDED];
        };
        uint32_t poly_word[ADDITIONAL+R_SIZE_PRIME+PADDED];
    };
}poly_buffer_t;

//ARM NEON Custom Cyclic shift instruction based binary polynomial shift left
void General_cyclic_shift_left_cNeon(const uint32_t a_in[ADDITIONAL+R_SIZE_PRIME+PADDED],uint16_t current_shift,uint32_t res[R_SIZE_PRIME+PADDED]){
    size_t vl=vsetvl_e32m1(4);//vl should be 4

    uint32_t* additional_chunk_addr=a_in;
    uint32_t* inspect_chunk_addr=a_in+4;

    uint16_t total_chunk_num=(R_SIZE_PRIME+PADDED)>>2;//each chunk contains 128 bit, for 12323bit there are 97 128-bit chunks
    uint8_t chunk_shift = current_shift >> 7;

    uint32_t* res_chunk_addr=res+(chunk_shift<<2);

    uint8_t bit_shift1 = current_shift & 127;
    uint8_t bit_shift2 = bit_shift1 + 128 - (R_BITS & 127);
    bool flag=bit_shift2>=128;
    bit_shift2=flag?bit_shift2-128:bit_shift2;

    uint16_t temp;//has no meaningful usage
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(temp):[s]"r"(bit_shift1),[imm]"i"(0));

    for(int i=0;i<total_chunk_num-chunk_shift;i++){
        vuint32m1_t v_additional=vle32_v_u32m1(additional_chunk_addr,vl);
        vuint32m1_t v_inspect=vle32_v_u32m1(inspect_chunk_addr,vl);
        vuint32m1_t v_cycres;
        __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(v_cycres):[vt]"vr"(v_inspect),[vs]"vr"(v_additional));
        vuint32m1_t v_res=vle32_v_u32m1(res_chunk_addr,vl);
        v_res=vxor_vv_u32m1(v_res,v_cycres,vl);
        vse32_v_u32m1(res_chunk_addr,v_res,vl);

        //addres update
        additional_chunk_addr+=4;
        inspect_chunk_addr+=4;
        res_chunk_addr+=4;
    }

    res_chunk_addr=res;
    if(flag){
        additional_chunk_addr-=4;
        inspect_chunk_addr-=4;
    }
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(temp):[s]"r"(bit_shift2),[imm]"i"(0));
    for(int i=0;i<chunk_shift;i++){
        vuint32m1_t v_additional=vle32_v_u32m1(additional_chunk_addr,vl);
        vuint32m1_t v_inspect=vle32_v_u32m1(inspect_chunk_addr,vl);
        vuint32m1_t v_cycres;
        __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(v_cycres):[vt]"vr"(v_inspect),[vs]"vr"(v_additional));
        vuint32m1_t v_res=vle32_v_u32m1(res_chunk_addr,vl);
        v_res=vxor_vv_u32m1(v_res,v_cycres,vl);
        vse32_v_u32m1(res_chunk_addr,v_res,vl);

        //addres update
        additional_chunk_addr+=4;
        inspect_chunk_addr+=4;
        res_chunk_addr+=4;
    }
}

//Only Simulate Arm NEON work procedure
int main(){
    poly_buffer_t dense_poly;
    uint16_t h1_compact[DV];
    for(int i=0;i<DV;i++)h1_compact[i]=rand()%R_BITS;

    uint32_t* cycsll_res=(uint32_t*)calloc(R_SIZE_PRIME+PADDED,sizeof(uint32_t));

    size_t vl=vsetvl_e32m1(4);//vl should be 4
    
    //Get Additional Chunk Process
    vuint32m1_t vx=vle32_v_u32m1(&dense_poly.poly_word[388],vl);
    vuint32m1_t vy=vle32_v_u32m1(&dense_poly.poly_word[384],vl);
    uint16_t addchunk_shift_offset=128-R_BITS & 127;
    uint16_t shift_offset;
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(shift_offset):[s]"r"(addchunk_shift_offset),[imm]"i"(0));
    __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(vy):[vt]"vr"(vx),[vs]"vr"(vy));
    vse32_v_u32m1(dense_poly.poly_word,vy,vl);

    //multiplication process
    unsigned long start, end;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    for(int i=0;i<DV;i++) General_cyclic_shift_left_cNeon(dense_poly.poly_word,h1_compact[i],cycsll_res);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("total time spent by poly_mul_RVV:%lu\n",end-start);

    //final mask
    cycsll_res[R_SIZE_PRIME-1]&=FINAL_CHUNK_MASK;
    cycsll_res[R_SIZE_PRIME]=cycsll_res[R_SIZE_PRIME+1]=0;

    //print out
    //for(int i=0;i<R_SIZE_PRIME;i++)printf("cyclic_shift_res[%d]=%u\n",i,cycsll_res[i]);

    free(cycsll_res);
}