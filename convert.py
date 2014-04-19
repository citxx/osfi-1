#!/usr/bin/env python3

f = open('absorption.csv')
values = [line.split(',') for line in f.readlines()]
f.close()

print(values)
mx = max(float(value) for freq, value in values)
new_values = [[freq, str(float(value) / mx)] for freq, value in values]

f = open('absorption2.csv', 'w')

for freq, value in new_values:
  print(freq, value, sep=',', file=f)

f.close()
