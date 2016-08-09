    x_depart = x - 0.5*dx - u*dt
    while len(np.where(np.logical_or(x_depart < x[0], x_depart > x[-1]))[0]) > 0:
        x_depart = np.where(x_depart < x[0], x_depart + L, x_depart)
        x_depart = np.where(x_depart > x[-1], x_depart - L, x_depart)
    # idx = np.where(idx ==0, -1, idx)
    x_r = np.zeros_like(x)
    x_r[:-1] = x[:-1] + 0.5*(x[1:] - x[:-1])
    x_r[-1] = x[-1] + 0.5*(x[1] - x[0])
    x_l = np.zeros_like(x) 
    x_l[1:] = x[1:] - 0.5*(x[1:] - x[:-1])
    x_l[0] = x[0] - 0.5*(x[-1] - x[-2])
    # idx = np.where(idx ==0, -1, idx)
    for i in xrange(nx+1):
        for j in xrange(nx+1):
            if  x_r[i] > x_depart[j] and  x_depart[j] > x_l[i]:
                idx[j] = i
    for i in xrange(nx+1):
        for j in xrange(nx+1):
            if abs(x_depart[j] - x_r[i])<10e-10:
                idx[j] = i