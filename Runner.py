from Schemes import *
from matplotlib import animation

schemes = (LF, LW, LF_cons, LW_cons, Godunov)
names = ["LF", "LW", "LF cons", "LW cons", "Godunov"]


p_vec = [np.copy(p0) for p in range(len(schemes))]
x = np.arange(0, 1, dx) 
p_vec.append(x)

fig = plt.figure()
ax1 = plt.axes(xlim=(0, 1), ylim=(0, .4))

lines = []

for index in range(len(schemes) + 1):
    if index == len(schemes):
        lobj = ax1.plot([],[], "ro", markersize=1, label="Lagrange")[0]
    else:
        lobj = ax1.plot([],[], label=names[index])[0]
    lines.append(lobj)


def init():
    plt.title("Comparison of different schemes")
    plt.xlabel("x")
    plt.ylabel("Rho")
    plt.legend()
    for line in lines:
        line.set_data([],[])
    return lines


def animate(i):
    for n in range(len(schemes)):
        p_vec[n] = schemes[n](p_vec[n])
        lines[n].set_data(x, p_vec[n])
    p_vec[len(schemes)] = Lagrange(p_vec[len(schemes)])
    lines[len(schemes)].set_data(p_vec[len(schemes)], p0)
    return lines

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1, interval=1, blit=True)

plt.show()




