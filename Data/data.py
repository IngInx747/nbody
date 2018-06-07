import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from random import randint

#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111)

plt.xlabel('Nodes')
#plt.ylabel('Runtime (s)')
plt.ylabel('Speedup')
t1=3.688e+02
t2=3.571e+02
x1=[32,64,128,256,512,1024]
y1=[
3.688e+02,
1.855e+02,
9.441e+01,
5.008e+01,
2.422e+01,
1.315e+01,]
y2=[
3.571e+02,
1.813e+02,
9.538e+01,
5.380e+01,
4.148e+01,
5.768e+01,]
y0 = [t/32 for t in x1]
y1 = [t1/y1[i] for i in range(len(y1))]
y2 = [t2/y2[i] for i in range(len(y2))]
plt.plot(x1,y0,'k-')
plt.plot(x1,y1,'bo-.',label="Hybrid")
plt.plot(x1,y2,'rs--',label="Pure")

leg = plt.legend(ncol=1)
leg.get_frame().set_alpha(0.5)

#ax.plot(points)
#plt.savefig('ss_6_BGQ.png')
plt.savefig('ss_6_BGQ_speedup.png')