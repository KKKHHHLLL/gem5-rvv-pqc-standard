#ifndef __ARCH_RISCV_INSTS_PQC_HH__
#define __ARCH_RISCV_INSTS_PQC_HH__

#include <string>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/regs/misc.hh"
#include "arch/riscv/utility.hh"
#include "arch/riscv/types.hh"
#include "cpu/exec_context.hh"
#include "cpu/static_inst.hh"

#include <assert.h>
#include "arch/riscv/insts/vector.hh"

namespace gem5
{

namespace RiscvISA
{

/**
 * Base class for Vector Config operations
 */

class PqcConfOp : public RiscvStaticInst
{
  protected:
    PqcConfOp(const char *mnem, ExtMachInst _extMachInst, OpClass __opClass)
        : RiscvStaticInst(mnem, _extMachInst, __opClass)
    {}

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

class PqcWholeVecArithOp : public RiscvStaticInst
{
  protected:
    PqcWholeVecArithOp(const char *mnem, ExtMachInst _extMachInst, OpClass __opClass)
        : RiscvStaticInst(mnem, _extMachInst, __opClass)
    {}

  std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

}

}

#endif // __ARCH_RISCV_REGS_PQC_HH__