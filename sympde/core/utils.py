# coding: utf-8

import re
import string
import random
from sympy.utilities.iterables import cartes

__all__ = ('random_string', 'expand_name_patterns')

#==============================================================================
def random_string( n ):
    chars    = string.ascii_uppercase + string.ascii_lowercase
    selector = random.SystemRandom()
    return ''.join( selector.choice( chars ) for _ in range( n ) )

#==============================================================================
_range = re.compile('([0-9]*:[0-9]+|[a-zA-Z]?:[a-zA-Z])')

def expand_name_patterns(names, *, seq=None):
    """
    Expand names according to the given pattern or list of patterns, following
    exactly the same rules of Sympy's "symbols" function.

    Parameters
    ----------
    names: str or iterable
        Variable name, or name pattern, or list of name patterns, or list of
        [list of ...] name patterns.

    seq: bool
        Set to True if an iterable container is needed for a single name.

    Results
    -------
    res: str or iterable
        Single name string, or list of [list of ...] name strings (having the
        same type as the input parameter 'names').

    See also
    --------
    https://docs.sympy.org/latest/modules/core.html#sympy.core.symbol.symbols

    Notes
    -----
    This code is extracted from Sympy's "symbols" function, but it only
    generates variables' names without creating new objects, which can be
    instantiated at a later stage. This creates a very clear logical
    separation between 1) names' generation, and 2) objects' instantiation,
    plus it allows one to reuse the same name expansion rules in different
    locations of the library, for generating different objects. Further,
    since separate functions take on conceptually different tasks, it is much
    easier to test them individually.

    This approach would also be beneficial to Sympy's "symbols", as it would
    allow for other Python libraries based on Sympy (like SymPDE) to be
    compliant with Sympy's name expansion rules without duplicating its code.
    Hence, in the future we would like to propose the integration of this
    function into Sympy.

    """
    result = []

    if isinstance(names, str):
        marker = 0
        literals = [r'\,', r'\:', r'\ ']
        for i in range(len(literals)):
            lit = literals.pop(0)
            if lit in names:
                while chr(marker) in names:
                    marker += 1
                lit_char = chr(marker)
                marker += 1
                names = names.replace(lit, lit_char)
                literals.append((lit_char, lit[1:]))
        def literal(s):
            if literals:
                for c, l in literals:
                    s = s.replace(c, l)
            return s

        names = names.strip()
        as_seq = names.endswith(',')
        if as_seq:
            names = names[:-1].rstrip()
        if not names:
            raise ValueError('no symbols given')

        # split on commas
        names = [n.strip() for n in names.split(',')]
        if not all(n for n in names):
            raise ValueError('missing symbol between commas')

        # split on spaces
        for i in range(len(names) - 1, -1, -1):
            names[i: i + 1] = names[i].split()

        # if provided, get 'seq' value from parameter list
        if seq is None:
            seq = as_seq
        else:
            if not isinstance(seq, bool):
                msg = 'type(seq) must be bool, got {} instead'.format(type(seq))
                raise TypeError(msg)

        # expand ranges
        for name in names:
            if not name:
                raise ValueError('missing symbol')

            if ':' not in name:
                result.append(literal(name))
                continue

            split = _range.split(name)
            # remove 1 layer of bounding parentheses around ranges
            for i in range(len(split) - 1):
                if i and ':' in split[i] and split[i] != ':' and \
                        split[i - 1].endswith('(') and \
                        split[i + 1].startswith(')'):
                    split[i - 1] = split[i - 1][:-1]
                    split[i + 1] = split[i + 1][1:]
            for i, s in enumerate(split):
                if ':' in s:
                    if s[-1].endswith(':'):
                        raise ValueError('missing end range')
                    a, b = s.split(':')
                    if b[-1] in string.digits:
                        a = 0 if not a else int(a)
                        b = int(b)
                        split[i] = [str(c) for c in range(a, b)]
                    else:
                        a = a or 'a'
                        split[i] = [string.ascii_letters[c] for c in range(
                            string.ascii_letters.index(a),
                            string.ascii_letters.index(b) + 1)]  # inclusive
                    if not split[i]:
                        break
                else:
                    split[i] = [s]
            else:
                seq = True
                if len(split) == 1:
                    names = split[0]
                else:
                    names = [''.join(s) for s in cartes(*split)]
                if literals:
                    result.extend([literal(s) for s in names])
                else:
                    result.extend([s for s in names])

        if not seq and len(result) <= 1:
            if not result:
                return ()
            return result[0]

        return tuple(result)

    else:
        # recursive call on list of names, seq is ignored
        for name in names:
            result.append(expand_name_patterns(name))

        return type(names)(result)
