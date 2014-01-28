
from abc import ABCMeta, abstractmethod, abstractproperty
import os
import time
import random
import math
import cStringIO


class AbstractGenerator(object):
  __metaclass__ = ABCMeta

  @abstractmethod
  def run(self, input_params):
    raise NotImplementedError( "Should have implemented this" )

  @abstractmethod
  def load(self, conn, d):
    raise NotImplementedError( "Should have implemented this" )

