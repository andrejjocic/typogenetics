"""Interpreter for typogenetics.
Binding site denoted by ↑ (or ↟ with copy mode enabled).

Valid instructions (apart from typo-acids):
- help: print this message
- halt: halt execution
- flip: toggle upside-down printing of complementary bases
- typ: print all strands in the typo-plasm
- sub: print all strands in the substrate
"""

from typogenetics import *


def debug_print():
    print(f"\nstart: {strand}@{init_unit}")
    print("applied enzyme:", enzyme) 


print(__doc__)
strand = Strand(input("enter initial strand: "))
init_unit = int(input(f"enter binding unit (1 to {len(strand)}): "))

print_flipped = True

env = EnzymeRuntime(substrate=Molecule.from_strand(strand, init_unit - 1),
                    interactive=True)
env.print()
enzyme: Enzyme = []

while True:
    try:
        enzyme.append(
            amino_acid := AminoAcid[instr := input("instruction: ").lower()]
        )
        try:
            env.execute(amino_acid)
            env.print(flip_upper=print_flipped)

        except EnzymeRuntime.Halt:
            print("the enzyme fell off.")
            break
        
        except:
            print("!!! bug in typogenetics.py !!!")
            debug_print()
            raise

    except KeyError: # non-aminoacid instructions
        match instr:
            case "help":
                print(__doc__)

            case "flip":
                print_flipped ^= True
                env.print(flip_upper=print_flipped)

            case "typ":
                print(sum((mol.get_strands() for mol in env.cytoplasm), start=[]))
 
            case "sub":
                strands = env.substrate.get_strands(active_idx=-1)
                print("active strand:", strands.pop())
                print("remaining strands:", strands)

            case "halt":
                break

            case _:
                print("invalid instruction:", instr)

    except KeyboardInterrupt:
        debug_print()
        break


try:
    strands = env.substrate.get_strands(active_idx=-1)
    print("last active strand:", strands.pop())

except NoActiveStrand as exc:
    print("the active strand was just deleted.")
    strands = exc.inactive_strands

print("inactive strands (including typo-plasm):",
      sum((mol.get_strands() for mol in env.cytoplasm), start=strands))