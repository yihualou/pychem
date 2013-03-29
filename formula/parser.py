from ply import lex, yacc
from models import Molecule, Formula, Equality

class FormulaParser(object):
  def parse(self, table, formula):
    tokens = ("SYMBOL", "COUNT", "ADD", "EQUALS", "LPAREN", "RPAREN")
    t_ADD = "\s*\+\s*"
    t_EQUALS = "\s*=\s*"
    t_LPAREN = "\s*\(\s*"
    t_RPAREN = "\s*\)\s*"

    @lex.TOKEN(self.__compute_symbols(table))
    def t_SYMBOL(t):
      t.value = table[t.value]
      return t

    def t_COUNT(t):
      r"\d+"
      t.value = int(t.value)
      return t

    def t_error(t):
      raise TypeError("Unknown text '%s'" % (t.value,))

    lex.lex()

    def p_equation(p):
      """
      chemical_equation : chemical_formula
      chemical_equation : chemical_formula EQUALS chemical_formula
      """
      if len(p) == 1:
        p[0] = p[1]
      else:
        p[0] = Equality(p[1], p[3])

    def p_chemical_formula(p):
      """
      chemical_formula : molecule
      chemical_formula : chemical_formula ADD molecule
      """
      if len(p) == 2:
        p[0] = Formula([p[1]])
      else:
        formula = p[1]
        formula.molecules.append(p[3])
        p[0] = formula

    def p_molecule_species(p):
      """
      molecule : COUNT LPAREN species_list RPAREN
      molecule : species_list
      """
      if len(p) == 5:
        p[0] = (Molecule(p[3]), p[1])
      else:
        p[0] = (Molecule(p[1]), 1)

    def p_species_list(p):
      """
      species_list : species_list species
      species_list : species
      """
      if len(p) == 3:
        p[0] = p[1] + [p[2]]
      else:
        p[0] = [p[1]]

    def p_species_enclosed(p):
      """
      species : LPAREN species_list RPAREN COUNT
      """
      p[0] = (Molecule(p[2]), p[4])

    def p_species(p):
      """
      species : SYMBOL
      species : SYMBOL COUNT
      """
      if len(p) == 2:
        p[0] = (p[1], 1)
      elif len(p) == 3:
        p[0] = (p[1], p[2])

    def p_error(p):
      print "Syntax error at '%s'" % p.value

    yacc.yacc()
    return yacc.parse(formula)

  def __compute_symbols(self, table):
    l = [abbr for abbr in table]
    l.sort(key=len, reverse=True)
    return "|".join(l)
