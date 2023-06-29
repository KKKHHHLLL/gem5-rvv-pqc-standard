#include "arch/riscv/insts/pqc.hh"
#include <sstream>
#include <string>
#include <algorithm>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/utility.hh"
#include "cpu/static_inst.hh"

#include "arch/riscv/regs/vector.hh"

namespace gem5
{

namespace RiscvISA
{

std::string
PqcConfOp::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    std::string left_bracket="(";
    std::string right_bracket=")";
    auto source0=registerName(srcRegIdx(0))+left_bracket+std::to_string((uint8_t)machInst.vs1)+right_bracket;
    auto source1=registerName(srcRegIdx(1))+left_bracket+std::to_string((uint8_t)machInst.vs2)+right_bracket;
    ss << mnemonic << ' ' << source0 << ", " << source1;
    return ss.str();
}

std::string
PqcWholeVecArithOp::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", "
        << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    return ss.str();
}

}

}