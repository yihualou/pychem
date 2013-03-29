import numpy as np
import itertools

from collections import defaultdict
from pychem.util.matrix import nullspace

class Molecule(object):
  def __init__(self, components):
    condensed = defaultdict(int)
    for el, count in components:
      condensed[el] += count

    # Tuple of element, count or molecule, count tuples.
    self.components = tuple((el, count) for el, count in condensed.iteritems())

  def __hash__(self):
    return self.components.__hash__()

  def __eq__(self, other):
    return frozenset(self.components) == frozenset(other.components)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __repr__(self):
    return "molecule:" + repr(self.components)

  def __len__(self):
    return self.components.__len__()

  def __iter__(self):
    return self.components.__iter__()

  def __getitem__(self, key):
    return self.components.__getitem__(key)

  def molecular_mass(self):
    s = 0
    for item, count in self.components:
      if type(item) == Molecule:
        s += item.molecular_mass() * count
      else:
        s += item.atomic_mass * count
    return s

  def element_stats(self, mult=1):
    s = defaultdict(int)
    for item, count in self.components:
      if type(item) == Molecule:
        for se, sc in item.element_stats(count).iteritems():
          s[se] += sc
      else:
        s[item] += count

    for k, v in s.iteritems():
      s[k] *= mult

    return s

class Formula(object):
  def __init__(self, molecules):
    # Tuple of molecule, count tuples.
    self.molecules = tuple(molecules)

  def __hash__(self):
    return self.molecules.__hash__()

  def __eq__(self, other):
    return frozenset(self.molecules) == frozenset(other.molecules)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __repr__(self):
    return "formula:" + repr(self.molecules)

  def __len__(self):
    return self.molecules.__len__()

  def __iter__(self):
    return self.molecules.__iter__()

  def __getitem__(self, key):
    return self.molecules.__getitem__(key)

class Equality(object):
  def __init__(self, left, right):
    # Mutable left and right formulas.
    self.left = left
    self.right = right

  def balance(self):
    lstats = [m.element_stats(c) for m, c in self.left]
    rstats = [m.element_stats(c) for m, c in self.right]
    offset = len(lstats)

    # Calculate element mapping
    mapping = {}
    i = 0
    for s in lstats:
      for e in s.iterkeys():
        if e not in mapping:
          mapping[e] = i
          i += 1

    # Initialize dimensions and matrix
    n = len(mapping)
    m = len(lstats) + len(rstats)

    mat = np.zeros((n, m))

    for j, s in enumerate(lstats):
      for e, c in s.iteritems():
        mat[mapping[e]][j] = c
        
    for j, s in enumerate(rstats):
      for e, c in s.iteritems():
        mat[mapping[e]][offset + j] = -c

    ns = nullspace(mat)
    l = ns.T[0].tolist()
    mns = min(l, key=lambda v: abs(v))

    coefs = self.__normalize([e / mns for e in l])
    
    nl = []
    nr = []

    for i in xrange(offset):
      m, oldc = self.left[i]
      nl.append((m, coefs[i]))

    for i in xrange(len(coefs) - offset):
      m, oldc = self.right[i]
      nr.append((m, coefs[offset + i]))

    self.left = Formula(nl)
    self.right = Formula(nr)

  def __normalize(self, l, tol=1e-6):
    done = False
    i = 0

    while not done:
      i += 1
      done = True

      for e in l:
        v = e * i
        if abs(round(v) - v) > tol:
          done = False

    return [int(round(e * i)) for e in l]

  def is_balanced(self):
    lcomb = defaultdict(int)
    rcomb = defaultdict(int)

    for m, c in self.left:
      for e, c2 in m.element_stats(c).iteritems():
        lcomb[e] += c2

    for m, c in self.right:
      for e, c2 in m.element_stats(c).iteritems():
        rcomb[e] += c2

    return lcomb == rcomb

   