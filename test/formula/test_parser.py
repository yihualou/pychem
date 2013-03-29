import unittest

from pychem.periodic.table import PeriodicTable
from pychem.formula.parser import FormulaParser
from pychem.formula.models import Formula, Equality, Molecule

class TestParser(unittest.TestCase):
  @classmethod
  def setUpClass(cls):    
    cls.pt = PeriodicTable("res/periodic_table.csv")
    cls.fp = FormulaParser()

  def test_parse(self):
    pt = TestParser.pt

    f = TestParser.fp.parse(pt, "CH3CH2CH2CH2CH2CH3")
    fe = Formula([
      (Molecule([(pt['H'], 14), (pt['C'], 6)]), 1)  
    ])
    self.assertEquals(f, fe)

    f = TestParser.fp.parse(pt, "Ag2SeO4")
    fe = Formula([
      (Molecule([(pt['Ag'], 2), (pt['Se'], 1), (pt['O'], 4)]), 1)  
    ])
    self.assertEquals(f, fe)

    f = TestParser.fp.parse(pt, "Ca3(PO4)2 + H2O")
    fe = Formula([
      (Molecule([
        (pt['Ca'], 3),
        (Molecule([
          (pt['P'], 1),
          (pt['O'], 4)
        ]), 2)
      ]), 1),
      (Molecule([
        (pt['H'], 2), (pt['O'], 1)
      ]), 1)  
    ])
    self.assertEquals(f, fe)
