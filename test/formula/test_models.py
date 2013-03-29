import unittest

from pychem.periodic.table import PeriodicTable
from pychem.formula.parser import FormulaParser
from pychem.util.matrix import nullspace

class TestEquality(unittest.TestCase):
  @classmethod
  def setUpClass(cls):    
    cls.pt = PeriodicTable("res/periodic_table.csv")
    cls.fp = FormulaParser()

  def test_balance(self):
    tests = [
      "C5H12 + O2 = CO2 + H2O",
      "Zn + HCl = ZnCl2 + H2",
      "Ca(OH)2 + H3PO4 = Ca3(PO4)2 + H2O",
      "FeCl3 + NH4OH = Fe(OH)3 + NH4Cl",
      "S8 + F2 = SF6",
      "C2H6 + O2 = CO2 + H2O",
      "Al2(CO3)3 + H3PO4 = AlPO4 + CO2 + H2O"
    ]

    coef = [
      [1, 8, 5, 6],
      [1, 2, 1, 1],
      [3, 2, 1, 6],
      [1, 3, 1, 3],
      [1, 24, 8],
      [2, 7, 4, 6],
      [1, 2, 2, 3, 3]
    ]

    for i, test in enumerate(tests):
      equality = TestEquality.fp.parse(TestEquality.pt, test)
      self.assertEquals(equality.balance(), coef[i])