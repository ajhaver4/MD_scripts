import sys
import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv(sys.argv[1],delimiter='\s+',names=['time','1','2','3','4','5','6','7'])

fig,ax = plt.subplots()

#ax.plot(data['time'],data['2'],label='d1,d2 < 12.5 nm')
ax.plot(data['time'],data['3'],label='d1 < 12.5 nm')
#ax.plot(data['time'],data['4'],label='d1,d2 < 16 nm')
ax.plot(data['time'],data['5'],label=' Cry ')
ax.plot(data['time'],data['6'],label='4 < d1,d2 < 8 nm')
ax.plot(data['time'],data['7'],label='d1<5 ; 6 < d2 < 12 nm')

ax.legend()
ax.set_xlabel(r'Time ($\mu$s)')
ax.set_ylabel(r'$\Delta G$ (kJ/mol)')
plt.show()
