#! /usr/bin/env python3

import sys
from collections import Counter

line_number = 0
content_counter = Counter()
for line in sys.stdin:
	line_number += 1
	if line_number % 4 == 2:
		base_counter = Counter(line)
		try:
			content_counter[round((
				(base_counter['G'] + base_counter['C']) /
				(base_counter['A'] + base_counter['C'] + base_counter['G'] + base_counter['T'])
			), 2)] += 1 # round to 2 digits
		except ZeroDivisionError:
			pass

for content in sorted(content_counter.keys()):
	print('%.2f\t%i' % (content, content_counter[content])) # print 2 digits

