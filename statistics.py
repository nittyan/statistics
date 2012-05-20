# -*- coding:utf-8 -*-

import math
import sys


t_distribution_diagram = {
  1: [12.706, 63.657],
  2: [4.3027, 9.9248],
  3: [3.1825, 5.8409],
  4: [2.7764, 4.6041],
  5: [2.5706, 4.0321],
  6: [2.4469, 3.7074],
  7: [2.3646, 3.4995],
  8: [2.3060, 3.3554],
  9: [2.2622, 3.2498],
  10: [2.2281, 3.1693],
  11: [2.2010, 3.1058],
  12: [2.1788, 3.0545],
  13: [2.1604, 3.0123],
  14: [2.1448, 2.9768],
  15: [2.1315, 2.9467],
  16: [2.1199, 2.9208],
  17: [2.1098, 2.8982],
  18: [2.1009, 2.8784],
  19: [2.0930, 2.8609],
  20: [2.0860, 2.8453],
  21: [2.0796, 2.8314],
  22: [2.0739, 2.8188],
  23: [2.0687, 2.8073],
  24: [2.0639, 2.7969],
  25: [2.0595, 2.7874],
  26: [2.0555, 2.7787],
  27: [2.0518, 2.7707],
  28: [2.0484, 2.7633],
  29: [2.0452, 2.7564],
  30: [2.0423, 2.7500],
  40: [2.0211, 2.7045],
  60: [2.0003, 2.6603],
  120: [1.9799, 2.6174],
  181: [1.9600, 2.5758]
}


def average(seq):
  return sum(seq) / float(len(seq))


def dispersion(seq):
  '''分散'''
  
  avg = average(seq)  
  total = sum(_diff_from_avg(seq, avg))
  
  return total / float(len(seq))


def standard_variation(seq):
  '''標準偏差'''

  return math.sqrt(dispersion(seq))


def unbiased_estimate_of_variance(seq):
  '''不偏分散'''
  
  avg = average(seq)
  total = sum(_diff_from_avg(seq, avg))
  
  return total / float(len(seq) - 1)


def _diff_from_avg(seq, avg):
  return [math.pow(math.fabs(n - avg), 2) for n in seq]


def confidence_interval(seq):
  '''信頼区間'''
  
  freedom = len(seq) - 1
  
  # とりあえず信頼区間95%で計算
  t = _select_t(freedom)
  min = 0
  max = 0
  
  se = math.sqrt(unbiased_estimate_of_variance(seq) / len(seq))
  
  print(se)
  print(t)
  min = average(seq) + t * (- 1) * se
  max = average(seq) + t * se
  
  return (min, max)
  

# interval 95とか
def _select_t(freedom, interval = 0):
  '''自由度に当てはまるtの値を取得します'''
  
  if t_distribution_diagram.has_key(freedom):
    return t_distribution_diagram[freedom][0]
  
  min_t = 0
  min = sys.maxint
  for t in t_distribution_diagram.keys():
    if math.fabs(freedom - t) < min:
      min = math.fabs(freedom - t)
      min_t = t

  return t_distribution_diagram[min_t][0]


def chi_square(actual, expect):
  '''カイ2乗値'''
  
  total = 0
  for i in range(len(actual)):
    row = actual[i]
    
    for j in range(len(row)):
      total += math.pow(math.fabs(row[j] - expect[i][j]), 2) / expect[i][j] 

  return total


def covariance(seq1, seq2):
  '''共分散'''
  if len(seq1) != len(seq2):
    raise InvalidArgumentError()

  avg1 = average(seq1)
  avg2 = average(seq2)
  
  total = 0
  for i in range(len(seq1)):
    total += (seq1[i] - avg1) * (seq2[i] - avg2)

  return total / float(len(seq1))


class InvalidArgumentError(Exception):
  pass

class TTest:
  
  @classmethod
  def diff_of_standard_error(cls, seq1, seq2):
    '''差の標準誤差'''
    return math.sqrt(TTest.illation_of_population(seq1, seq2) * (1.0 / len(seq1) + 1.0 / len(seq2)))
    
  @classmethod
  def illation_of_population(cls, seq1, seq2):
    '''推定母分散'''

    return (dispersion(seq1) * len(seq1) + dispersion(seq2) * len(seq2)) / ((len(seq1) - 1) + (len(seq2) - 1))
  
  @classmethod
  def confidence_interval(cls, seq1, seq2): 
    '''信頼区間'''
    pass
    return average(seq1) - average(seq2) - 2


    
if __name__ == '__main__':
  wakuwaku = [
  3.5, 4.2, 4.9, 4.6, 2.8, 5.6, 4.2, 4.9, 4.4, 3.7,
  3.8, 4.0, 5.2, 3.9, 5.6, 5.3, 5.0, 4.7, 4.0, 3.1,
  5.8, 3.6, 6.0, 4.2, 5.7, 3.9, 4.7, 5.3, 5.5, 4.7,
  6.4, 3.8, 3.9, 4.2, 5.1, 5.1, 4.1, 3.6, 4.2, 5.0,
  4.2, 5.2, 5.3, 6.4, 4.4, 3.6, 3.7, 4.2, 4.8
  ]
  print(average(wakuwaku))
  print(dispersion(wakuwaku))
  print(standard_variation(wakuwaku))
  
  print('---------------------------')

  sample2 = [47, 51, 49, 50, 49, 46, 51, 48, 52, 49]
  print(unbiased_estimate_of_variance(sample2))
  
  
  print(confidence_interval(sample2))  
  print(chi_square([[435, 165],[265, 135]], [[420, 180],[280, 120]]))
  
  
  waku = [70, 75, 70, 85, 90, 70, 80, 75]
  mogu = [85, 80, 95, 70, 80, 75, 80, 90]
  
  print("illation")
  print(TTest.diff_of_standard_error(waku, mogu))

  print("共分散")
  print(covariance([71, 34, 58, 41, 69, 64, 16, 59], [64, 48, 59, 51, 56, 65, 45, 59]))

