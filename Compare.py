from Schemes import *

schemes = [LF, LW, LF_cons, LW_cons, Godunov]
names = ["LF", "LW", "LF cons", "LW cons", "Godunov"]

p_vec = [np.copy(p0) for p in range(len(schemes))]
x = np.arange(0, 1, dx) 
p_vec.append(x)

def L2(A, B):
    return np.sum(np.sqrt(((A - B)/A)**2))

def L_inf(A, B):
    return np.max((np.sqrt(((A - B)/A)**2)))

def check_collapse(x):
    test = (x - np.roll(x, -1)) > 0
    count= 0
    for value in test:
        if value:
            count += 1
        if count == 2:
            return True
    return False
    

time = []
def compare(T, norm=L2):
    t_collapse = 0
    n = len(schemes)
    differences = [[] for i in range(n)]


    for t in np.arange(0, T, dt):
        p_vec[n] = Lagrange(p_vec[n])
        lagrange_p = np.interp(x, p_vec[-1], p0, period=1)



        for i in range(n):
            p_vec[i] = schemes[i](p_vec[i])
            differences[i].append(norm(lagrange_p, p_vec[i]))
        time.append(t)
        if check_collapse(p_vec[-1]) and t_collapse == 0:
            #break
            t_collapse = t
        
    plt.title("Relative L2 norm")
    for i in range(n):
        plt.plot(time, differences[i], label=names[i])
    if t_collapse != 0:
        plt.plot([t_collapse, t_collapse], [0, np.max(differences)], '--')
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("Relative error")
    plt.show()

compare(.9, L_inf)
