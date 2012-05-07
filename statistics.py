# -*- coding:utf-8 -*-

import math


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
  190: [1.9600, 2.5758]
}


def average(seq):
  return sum(seq) / float(len(seq))


def dispersion(seq):
  '''分散'''
  
  avg = average(seq)  
  total = sum(diff_from_avg(seq, avg))
  
  return total / float(len(seq))


def standard_variation(seq):
  '''標準偏差'''

  return math.sqrt(dispersion(seq))


def unbiased_estimate_of_variance(seq):
  '''不偏分散'''
  avg = average(seq)
  total = sum(diff_from_avg(seq, avg))
  
  return total / float(len(seq) - 1)


def diff_from_avg(seq, avg):
  return [math.pow(math.fabs(n - avg), 2) for n in seq]


def confidence_interval():
  '''信頼区間'''



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
  