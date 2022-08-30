"""classes for simulating typogenetics"""

from __future__ import annotations
from typing import *
from enum import Enum, IntEnum, auto
from dataclasses import dataclass
import itertools
import math


AminoAcid = Enum("AminoAcid", start=0,
    # SEP: gene separator
    names="""
    SEP cut del swi
    mvr mvl cop off
    ina inc ing int
    rpy rpu lpy lpu
    """ # NOTE: names are in (row-major) order of typogenetic code
)
AminoAcid.__repr__ = lambda self: self.name

Enzyme = List[AminoAcid]


class Base(IntEnum):
    A = 0
    C = auto()
    G = auto()
    T = auto()

    @property
    def complement(self) -> Base:
        match self.name:
            case "A": return Base.T
            case "T": return Base.A
            case "G": return Base.C
            case "C": return Base.G

    @property
    def purine(self) -> bool:
        match self.name:
            case "A" | "G": return True
            case _: return False

    @property
    def pyrimidine(self) -> bool:
        return not self.purine

    def __repr__(self) -> str:
        return self.name

    __str__ = __repr__


@dataclass
class NoActiveStrand(Exception):
    """a function returning a list of strands was supposed to put
    the active one at a particular index, but there is no such strand"""
    
    inactive_strands: List[Strand]
    """the molecule's strands"""


class Strand:
    """DNA strand"""

    bases: List[Base]

    __hash: int

    def __init__(self, bases: List[Base] | str = None):
        if bases is None:
            self.bases = []

        elif isinstance(bases, str):
            cleaned = "".join(bases.split()).upper()
            try:
                self.bases = [Base[c] for c in cleaned]
            except KeyError:
                raise ValueError("invalid strand: " + cleaned)
        else:
            self.bases = bases

        self.compute_hash()


    def compute_hash(self) -> None:
        self.__hash = 0
        for b in self.bases: self.update_hash(b)


    BITS_PER_BASE = 3
    MAX_64BIT_LEN = math.ceil(64 / BITS_PER_BASE) - 1
    """max strand length such that its hash fits into 64 bits"""
    # not so important for Python, useful for
    # potential translation into another language

    def update_hash(self, last_base: Base) -> None:
        self.__hash <<= Strand.BITS_PER_BASE
        self.__hash |= last_base + 1
        
            
    def __eq__(self, obj):
        return hash(self) == hash(obj) and isinstance(obj, Strand)


    @property
    def translation(self) -> Iterator[Enzyme]:
        """yields the enzymes encoded in this strand"""

        curr_enzyme = []
        
        for i in range(0, len(self) - 1, 2):
            fst, snd = self[i:i + 2]

            if (acid := AminoAcid(4 * fst + snd)) is AminoAcid.SEP:
                if curr_enzyme:
                    yield curr_enzyme
                    curr_enzyme = []
            else:
                curr_enzyme.append(acid)
        
        if curr_enzyme: yield curr_enzyme


    class EnzymeNotApplicable(Exception): pass


    def apply_rightmost(self, enzyme: Enzyme) -> List[Strand]:
        """Returns the produced strands, with the active one at the start of the list.
        Assumes deterministic binding to rightmost unit with preferred base."""

        pref_base = binding_preference(enzyme)

        for i in range(len(self) - 1, -1, -1):
            if self[i] == pref_base: break
        else:
            raise Strand.EnzymeNotApplicable()

        runtime = EnzymeRuntime(
            substrate=Molecule.from_strand(self, init_pos=i)
        )
        for amino_acid in enzyme:
            try: runtime.execute(amino_acid)
            except EnzymeRuntime.Halt: break

        try:
            substrate_strands = runtime.substrate.get_strands(active_idx=0)
            return sum((mol.get_strands() for mol in runtime.cytoplasm),
                       start=substrate_strands)

        except NoActiveStrand as exc:
            for mol in runtime.cytoplasm:
                exc.inactive_strands += mol.get_strands()
            raise


    def apply(self, enzyme: Enzyme) -> Iterator[Tuple[List[Strand], bool]]:
        """Yields successive results of enzyme application, one for each initial binding.
        A result is a tuple (list of produced strands, flag) where the flag is True iff
        the active strand is at the start of the list (otherwise it isn't present, as it
        was deleted by the enzyme)."""
        
        pref_base = binding_preference(enzyme)

        for i, base in enumerate(self):
            if base == pref_base:
                runtime = EnzymeRuntime(
                    substrate=Molecule.from_strand(self, init_pos=i)
                )
                for amino_acid in enzyme:
                    try: runtime.execute(amino_acid)
                    except EnzymeRuntime.Halt: break

                try:
                    substrate_strands = runtime.substrate.get_strands(active_idx=0)
                    strands = sum((mol.get_strands() for mol in runtime.cytoplasm),
                                  start=substrate_strands)
                    yield strands, True
                    
                except NoActiveStrand as exc:
                    strands = sum((mol.get_strands() for mol in runtime.cytoplasm),
                                  start=exc.inactive_strands)
                    yield strands, False


    @property
    def daughters(self) -> List[Strand]:
        """Retruns the list (multiset) of daugher strands, resulting from self-application.
        Assumes enzyme binding to rightmost unit with preferred base."""

        dtrs = []
        active_strand = self

        for enzyme in self.translation:
            try:
                # OPT last call needn't place active strand at start
                strands = active_strand.apply_rightmost(enzyme)
                active_strand = strands[0]
                dtrs += strands[1:]
                    
            except NoActiveStrand as exc:
                dtrs += exc.inactive_strands
                return dtrs

            except Strand.EnzymeNotApplicable:
                dtrs.append(active_strand)
                return dtrs

        dtrs.append(active_strand)
        return dtrs


    @property
    def daughter_sets(self) -> Iterator[Iterator[Strand]]:
        """Yields successive results (strand multisets) of self-application,
        one for each sequence of initial enzyme bindings.
        
        WARNING: a yielded generator becomes invalid after the next one is yielded!"""

        enzymes = list(self.translation)
        
        # intermediate enzyme application results
        strand_groups: List[List[Strand]] = [None] * len(enzymes)

        def apply_enzymes_from(enzyme_idx: int, active_strand: Strand) \
                               -> Iterator[Iterator[Strand]]:

            if enzyme_idx == len(enzymes):
                yield itertools.chain([active_strand], *strand_groups)
            else:
                app_results = active_strand.apply(enzymes[enzyme_idx])
                try:
                    first = next(app_results)
                    for strands, active_present in itertools.chain([first], app_results):
                        if active_present:
                            new_active = strands[0]
                            strand_groups[enzyme_idx] = strands[1:]
                            yield from apply_enzymes_from(enzyme_idx + 1, new_active)
                        else:
                            yield itertools.chain(strands, *strand_groups[:enzyme_idx])

                except StopIteration: # enzyme not applicable to active strand
                    yield itertools.chain([active_strand], *strand_groups[:enzyme_idx])


        yield from apply_enzymes_from(0, active_strand=self)
        

    def reverse(self) -> None:
        self.bases.reverse()
        self.compute_hash()

    def __iadd__(self, base: Base):
        self.bases.append(base)
        self.update_hash(last_base=base)
        return self

    def __repr__(self):
        return "".join(b.name for b in self.bases) if self else 'ε'

    def __lt__(self, obj):
        return True # arbitrary value for tie-breaking in heap

    def __getitem__(self, key):
        return self.bases[key]

    def __hash__(self): return self.__hash
    def __len__ (self): return len(self.bases)
    def __bool__(self): return bool(self.bases)
    def __iter__(self): return iter(self.bases)



def binding_preference(enzyme: Enzyme) -> Base: 
    class Direction(IntEnum): # NOTE: clockwise
        RIGHT = 0
        DOWN = auto()
        LEFT = auto()
        UP = auto()

    last_link_direction = Direction.RIGHT
    
    for amino_acid in enzyme[1:-1]:
        match amino_acid.name:
            case "off" | "int" | "rpu" | "lpy" | "lpu":
                # turn left (couter-clockwise)
                last_link_direction = Direction((last_link_direction - 1) % len(Direction))

            case "swi" | "cop" | "inc" | "ing" | "rpy":
                # turn right (clockwise)
                last_link_direction = Direction((last_link_direction + 1) % len(Direction))
        
    match last_link_direction:
        case Direction.RIGHT: return Base.A
        case Direction.UP   : return Base.C
        case Direction.DOWN : return Base.G
        case Direction.LEFT : return Base.T



class Pair:
    """DNA molecule position, composed of 2 complementary units.
    At least one is always occupied by a base"""

    molecule: Molecule
    """the DNA molecule this unit is a part of"""

    # these names make sense if molecule isn't fliped (upside down)
    __lower_base: Base | None = None
    __upper_base: Base | None = None
    __left_pair: Pair | None = None
    __right_pair: Pair | None = None


    def __init__(self, containing_molecule: Molecule):
        self.molecule = containing_molecule


    @property
    def lower_base(self) -> Base | None:
        return self.__upper_base if self.molecule.flipped else self.__lower_base

    @lower_base.setter
    def lower_base(self, new_lower: Base | None) -> None:
        if self.molecule.flipped:
            self.__upper_base = new_lower
        else:
            self.__lower_base = new_lower


    @property
    def upper_base(self) -> Base | None:
        return self.__lower_base if self.molecule.flipped else self.__upper_base

    @upper_base.setter
    def upper_base(self, new_upper: Base | None) -> None:
        if self.molecule.flipped:
            self.__lower_base = new_upper
        else:
            self.__upper_base = new_upper


    def __repr__(self) -> str:
        return f"(low={self.lower_base}, upp={self.upper_base})"


    @property
    def left(self) -> Pair:
        """the left neighbor in the containing molecule"""
        return self.__right_pair if self.molecule.flipped else self.__left_pair

    @left.setter
    def left(self, new_left: Pair) -> None:
        if self.molecule.flipped:
            self.__right_pair = new_left
        else:
            self.__left_pair = new_left


    @property
    def right(self) -> Pair:
        """the right neighbor in the containing molecule"""
        return self.__left_pair if self.molecule.flipped else self.__right_pair

    @right.setter
    def right(self, new_right: Pair) -> None:
        if self.molecule.flipped:
            self.__left_pair = new_right
        else:
            self.__right_pair = new_right


    def disconnected_from(self, pair: Pair | None) -> bool:
        """assuming these are neighboring pairs, do they need to be unlinked?"""
        if pair is None: return False

        upper_disc = self.upper_base is None or pair.upper_base is None
        lower_disc = self.lower_base is None or pair.lower_base is None
        return upper_disc and lower_disc


    def unlink_left(self) -> Pair:
        """returns the disconnected neighbor"""

        if (old_left := self.left) is not None:
            self.left.right = None
            self.left = None

        return old_left


    def unlink_right(self) -> Pair:
        """returns the disconnected neighbor"""

        if (old_right := self.right) is not None:
            self.right.left = None
            self.right = None

        return old_right



@dataclass
class Molecule:
    """DNA molecule, consisting of two complementary (sparse) DNA sequences.
    It's either being acted on by an enzyme, or floating in the cytoplasm."""

    binding_site: Pair
    """the pair where an enzyme is bound (or an arbitrary pair, if there is no enzyme)"""

    flipped: bool = False


    @staticmethod
    def from_strand(strand: Strand, init_pos: int) -> Molecule:
        if len(strand) == 0:
            raise ValueError("empty initial strand not handled")

        if not (0 <= init_pos < len(strand)):
            raise ValueError(
                f"invalid index for strand of length {len(strand)}: {init_pos}")

        prev_pair = None
        mol = Molecule(binding_site=None)

        for base in strand:
            new_pair = Pair(containing_molecule=mol)
            new_pair.lower_base = base

            if prev_pair is not None: prev_pair.right = new_pair
            new_pair.left = prev_pair

            prev_pair = new_pair

        
        mol.binding_site = new_pair # rightmost

        for _ in range(len(strand) - init_pos - 1):
            mol.binding_site = mol.binding_site.left

        return mol


    def flip(self) -> None:
        self.flipped = not self.flipped


    def get_strands(self, active_idx: int = None) -> List[Strand]:
        """All the connected DNA strands in this molecule.
        The active one is placed at active_idx (if given)."""

        def split_strands(upper: bool) -> Iterator[Strand]:
            nonlocal active_strand

            base_getter = (lambda pair: pair.upper_base) if upper \
                      else lambda pair: pair.lower_base

            pairs = self.pairs
            if upper: pairs = reversed(list(pairs))

            curr_strand = Strand()

            for pair in pairs:
                if (base := base_getter(pair)) is not None:
                    curr_strand += base
                elif curr_strand:
                    yield curr_strand
                    curr_strand = Strand()

                if pair is self.binding_site and curr_strand:
                    active_strand = curr_strand

            # handle remainder
            if curr_strand: yield curr_strand

        
        active_strand = None
        
        # NOTE: order of calls important!
        # (lower strands are preferred when looking for active one)
        upper = list(split_strands(upper=True))
        lower = list(split_strands(upper=False))
        strands = upper + lower

        if active_idx is not None:
            try:
                active_at = strands.index(active_strand)
                try:
                    strands[active_idx], strands[active_at] = \
                        strands[active_at], strands[active_idx]
                except IndexError:
                    print(f"invalid index for list of {len(strands)} strands: {active_idx}")

            except ValueError: # strands.index failed
                raise NoActiveStrand(strands)
            
        return strands

    
    @property
    def pairs(self) -> Iterator[Pair]:
        if (p := self.binding_site) is None: return

        while p.left is not None:
            p = p.left

        while p is not None:
            yield p
            p = p.right
    


class EnzymeRuntime:
    """the environment of an active enzyme"""

    substrate: Molecule
    """the molecule the enzyme is bound to"""

    copy_mode: bool 

    cytoplasm: List[Molecule]
    """the molecule fragments that fell off during operation"""

    print_fallen: bool
    """print strands as they fall into the cytoplasm"""

    def __init__(self, substrate: Molecule, interactive: bool = False):
        self.substrate = substrate
        self.copy_mode = False
        self.cytoplasm = []
        self.print_fallen = interactive


    class Halt(Exception):
        """The enzyme fell off the molecule it was bound to.
        The molecule's binding site is on the last strand where the enzyme was bound
        (or None, if it was just deleted)."""
        pass


    @property
    def location(self) -> Pair:
        return self.substrate.binding_site

    @location.setter
    def location(self, new_loc: Pair) -> None:
        self.substrate.binding_site = new_loc

    
    def cond_copy(self, searching: bool) -> bool:
        """Account for copy mode (if enabled) on the current location.
        Returns True iff the lower base was filled in."""
        
        if self.copy_mode:
            if searching and self.location.lower_base is None:
                self.location.lower_base = self.location.upper_base.complement
                return True
            else:
                self.location.upper_base = self.location.lower_base.complement

        return False


    def discard(self, unlinked_pair: Pair) -> None:
        self.cytoplasm.append(mol := Molecule(binding_site=unlinked_pair))

        if self.print_fallen and (loose := mol.get_strands()):
            print("strands fell into typo-plasm:", loose)


    def move_left(self, searching: bool = False) -> bool:
        """return True iff the lower base was filled in"""

        if (left := self.location.left) is None or \
            not searching and left.lower_base is None:
            raise EnzymeRuntime.Halt()
        else:
            self.location = left
            return self.cond_copy(searching)


    def move_right(self, searching: bool = False) -> bool:
        """return True iff the lower base was filled in"""

        if (right := self.location.right) is None or \
            not searching and right.lower_base is None:
            raise EnzymeRuntime.Halt()
        else:
            self.location = right
            return self.cond_copy(searching)


    def cut_right(self) -> None:
        if (old_right := self.location.right) is None: return

        self.location.unlink_right()
        self.discard(old_right)


    def delete_base(self) -> None:
        """delete current base and move right"""

        self.location.lower_base = None
        old_left = self.location.left

        try:
            self.move_right()
            if self.location.disconnected_from(self.location.left):
                 self.discard(self.location.unlink_left())
            elif self.location.left.disconnected_from(old_left):
                self.discard(self.location.left.unlink_left())

        except EnzymeRuntime.Halt:
            # the last active strand is on the left of the deleted base
            # (if there's anything there)
            if (left := self.location.left) is not None and left.lower_base is not None:
                self.location = left
                if self.location.disconnected_from(self.location.right):
                    self.discard(self.location.unlink_right())
            else:
                self.discard(self.location)
                self.location = None
            raise

        


    def insert_right(self, base: Base) -> None:
        """Insert the base on the right and move to it"""

        new_pair = Pair(containing_molecule=self.substrate)
        new_pair.lower_base = base
        
        old_right = self.location.right
        self.location.right = new_pair
        new_pair.left = self.location

        self.move_right()
        
        if old_right is not None:
            self.location.right = old_right
            old_right.left = self.location

            if self.location.disconnected_from(old_right):
                self.discard(self.location.unlink_right())

        
    def execute(self, amino_acid: AminoAcid) -> None:
        match amino_acid.name:
            case "cop":
                self.copy_mode = True
                self.location.upper_base = self.location.lower_base.complement

            case "off":
                self.copy_mode = False

            case "swi":
                if self.location.upper_base is None:
                    raise EnzymeRuntime.Halt()
                else:
                    self.substrate.flip()

            case "cut":
                self.cut_right()

            case "del":
                self.delete_base()

            case insert if insert.startswith("in"):
                base = Base[insert[-1].upper()]
                self.insert_right(base)

            case move if insert.startswith("mv"):
                match move[-1]:
                    case 'l': self.move_left()
                    case 'r': self.move_right()

            case search:
                # OPT cache location of nearest pyri/puri (and just jump if copy=off)
                match search[0]:
                    case 'l': move_fun = self.move_left
                    case 'r': move_fun = self.move_right

                match search[1:]:
                    case "py": target_type = lambda base: base.pyrimidine
                    case "pu": target_type = lambda base: base.purine

                while True:
                    lower_filled = move_fun(searching=True)
                    
                    if (low_base := self.location.lower_base) is not None \
                            and not lower_filled and target_type(low_base):
                        break


    def print(self, flip_upper=True) -> None:
        pairs = list(self.substrate.pairs)

        strand_chars = lambda get_base: \
            [' ' if (b := get_base(pair)) is None else b.name for pair in pairs]

        upper_chars = strand_chars(lambda pair: pair.upper_base)
        lower_chars = strand_chars(lambda pair: pair.lower_base)
        
        bindsite_idx = pairs.index(self.location)

        if flip_upper:
            UPSIDE_DOWN = {'A': '∀', 'C': 'Ↄ', 'G': '⅁', 'T': '⊥'}
            upper_chars = [UPSIDE_DOWN[c] if c in UPSIDE_DOWN else ' ' for c in upper_chars]

        rows = [
            "_" * len(pairs),
            "".join(upper_chars),
            "".join(lower_chars),
            " " * bindsite_idx + ('↟' if self.copy_mode else '↑'),
            "‾" * len(pairs)
        ]
        for r in rows: print(r)
