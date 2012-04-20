from numpy import *
from sdf import *
import matplotlib.pyplot as plt

a=SDF("0010.sdf")
b=a.read()
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(b['grid/X'],b['ex'][:,64],'r+')
plt.show()

