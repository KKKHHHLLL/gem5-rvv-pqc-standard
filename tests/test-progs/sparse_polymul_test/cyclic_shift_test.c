#include <stdio.h>
#include <stdint.h>
#include <riscv_vector.h>

#define VLEN 128// in bit

//test for cycsll
int main() {
	uint8_t p1[VLEN >> 3];
	uint8_t p2[VLEN >> 3];
	uint8_t pd[VLEN >> 3];
	for (int i = 0; i < (VLEN >> 3); i++) {
		p1[i] = p2[i] = i;
	}
	uint16_t shift_offset = 45;
    uint32_t returned_rd;
	vuint8m1_t vx,vy,vz;
    size_t avl=VLEN>>3;
    size_t vl;
    vl=vsetvl_e8m1(avl);
    vx=vle8_v_u8m1(p1,vl);
    vy=vle8_v_u8m1(p2,vl);
    __asm__ __volatile__ ("cycshift_conf %[d],%[s],%[imm]":[d]"=r"(returned_rd):[s]"r"(shift_offset),[imm]"i"(0));
    printf("returned_rd=%d\n",returned_rd);
    __asm__ __volatile__ ("cycsll %[vd],%[vt],%[vs]":[vd]"=vr"(vz):[vt]"vr"(vy),[vs]"vr"(vx));
    vse8_v_u8m1(pd,vz,vl);
    for(int i=0;i<(VLEN>>3);i++)printf("pd[%d]=%d\n",i,pd[i]);
}