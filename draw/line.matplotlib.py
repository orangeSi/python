#!/usr/bin/env python3
## close xming
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("in.txt", sep='\t') #header is default, can change to header=None
plt.plot(data['A'])
plt.savefig('in.line.png')
