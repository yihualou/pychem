import csv

from pychem.util.singleton import Singleton
from element import Element

class PeriodicTable(object):
  __metaclass__ = Singleton

  def __init__(self, path):
    self.elements = {}

    with open(path, 'rb') as f:
      reader = csv.reader(f, delimiter=',', quotechar='|')
      for row in reader:
        if row:
          e = Element(row[2], row[3], int(row[0]), float(row[1]))
          self.elements[row[3]] = e

  def __len__(self):
    return self.elements.__len__()

  def __iter__(self):
    return self.elements.__iter__()

  def __getitem__(self, key):
    return self.elements.__getitem__(key)