import matplotlib.pyplot as plt
import pandas as pd

names = ['MCSteps', 'Energy', 'Average_Energy', 'Variance', 'Order_Parameter', 'Acceptance_Ratio', 'Temperature']
markers = ["r:", "k--", "b-", "g-."]
temp = [0.1, 0.5, 1]
# plot the LJ
fig1 = plt.figure(figsize =(14,7))
ax1 = fig1.add_subplot(111)
ax1.set_xlabel("Monte Carlo Steps", fontsize = 14)
ax1.set_ylabel(names[2], fontsize = 14)
ax1.set_title("LJ potentials as Function of Monte Carlo Steps",
              fontsize = 14,fontweight='bold')
ax1.grid(lw =0.5)

# plot the variance
fig2 = plt.figure(figsize =(14,7))
ax2 = fig2.add_subplot(111)

# plot the Order Parameter
fig3 = plt.figure(figsize =(14,7))
ax3 = fig3.add_subplot(111)

# plot the acceptance rate
fig4 = plt.figure(figsize =(14,7))
ax4 = fig4.add_subplot(111)

for i in range(3):
    df = pd.read_csv("bin/fcc{}.txt".format(i), names = names, header=0, sep = '\t')
    df.Order_Parameter = df.Order_Parameter.apply(lambda x : x if x < 1 else 1)
    if( i ==0):
         xp = ax1.twinx()
         xp.plot(df[names[2]],  markers[i], label ="T= {}".format(temp[i]), lw = 4)
         xp.set_ylabel(names[2] + " for T = {}".format(temp[i]), fontsize = 16, color = 'r')
         for tl in xp.get_yticklabels():
             tl.set_color('r')
             tl.set_fontsize(12)
    else:
        ax1.plot(df[names[2]],  markers[i], label ="T= {}".format(temp[i]), lw = 4)
        ax1.set_ylabel(names[2] + " for T = {} and T ={}".format(temp[1], temp[2]), fontsize = 16, color = 'k')
        for t in ax1.get_yticklabels():
             #t.set_color('r')
             t.set_fontsize(12)

    ax2.plot(df[names[3]],  markers[i], label ="T= {}".format(temp[i]), lw = 2, markersize = 1)
    ax3.plot(df[names[4]], markers[i], label ="T= {}".format(temp[i]), lw =2)
    ax4.plot(df[names[5]], markers[i], label ="T= {}".format(temp[i]), lw = 2)
    #print df.Average_Energy[0] #df.Order_Parameter[0]
    
ax1.legend(fontsize = 14, loc="center right", ncol=3)
#ax1.set_xlim(xmin = -250, xmax = 10250)

ax2.set_xlabel("Monte Carlo Steps", fontsize = 14)
ax2.set_ylabel(names[3], fontsize = 14)
ax2.set_title("Variances of LJ potentials as Function of Monte Carlo Steps",
              fontsize = 14,fontweight='bold')
ax2.legend(fontsize = 14, ncol=3)
#ax2.set_xlim(xmin = -250, xmax = 10250)
ax2.grid(lw = 0.5)


ax3.set_xlabel("Monte Carlo Steps", fontsize = 14)
ax3.set_ylabel(names[4], fontsize = 14)
ax3.set_title("Order Parameter as Function of Monte Carlo Steps",
              fontsize = 14,fontweight='bold')
ax3.legend(fontsize = 14, loc="lower left", ncol=3)
#ax3.set_xlim(xmin = -250, xmax = 10250)
ax3.set_ylim(ymax = 1.01)
ax3.grid(lw = 0.5)

ax4.set_xlabel("Monte Carlo Steps", fontsize = 14)
ax4.set_ylabel(names[5] +" %", fontsize = 14)
ax4.set_title("Acceptance Rate as Function of Monte Carlo Steps",
              fontsize = 14,fontweight='bold')
ax4.legend(fontsize = 14, loc="upper left", ncol=3)
#ax4.set_xlim(xmin = -250, xmax = 10250)
ax4.grid(lw = 0.5)

fig1.savefig("Potential.png")
fig2.savefig("Variance.png")
fig3.savefig("OrderParam.png")
fig4.savefig("AcceptanceRate.png")
plt.show()
