class Element(object):
  def __init__(self, full_name, abbr, atomic_number, atomic_mass):
    self.full_name = full_name
    self.abbr = abbr
    self.atomic_number = atomic_number
    self.atomic_mass = atomic_mass

  def __repr__(self):
    return self.abbr