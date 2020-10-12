#! /usr/bin/env python3

import sys
from collections import Counter

results = Counter()
line_no = 0
for line in sys.stdin:
	line_no += 1
	if line_no % 4 == 2:
		results[len(line.rstrip())] += 1

for length, count in sorted(results.items(), key = lambda item: item[0]): print(str(length) + '\t' + str(count))

