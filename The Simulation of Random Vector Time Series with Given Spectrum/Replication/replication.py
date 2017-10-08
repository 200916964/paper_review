import cmath
import numpy as np
# from matplotlib.pyplot import subplots
from numpy.linalg import inv, eigh

n = 256
q = 2
SIMULATION = 1000
OMEGA = np.array([[1.0, -0.5], [-0.5, 1.0]])
d = np.array([0.1, 0.4])
HALF_LAMBDA = 2.0 * np.pi / n * np.arange(start = 1, stop = int(n / 2) + 1, step = 1)
LAMBDA = 2.0 * np.pi / n * np.arange(start = 0, stop = n + 1, step = 1)


def D(L):
    diag = []
    for i in np.arange(0, q):
        diag.append(cmath.exp(d[i] * cmath.log(1 - L)))
    return np.diag(diag)


def F(lam):
    F = 0.5 * np.pi ** (-1) * np.matmul(np.matmul(inv(D(L = cmath.exp((-lam) * 1j))), OMEGA), inv(D(L = cmath.exp(lam * 1j))))
    # F[np.diag_indices_from(F)] = np.real(F[np.diag_indices_from(F)])
    return F


# def Star_Operation(v):
#     return np.transpose(np.conjugate(v))


def Randomise_Z(mean = 0, sd = 1):
    res = []
    re, im = np.random.normal(loc = mean, scale = sd, size = (2, 1)), np.random.normal(loc = mean, scale = sd, size = (2, 1))
    for i in np.arange(0, 2):
        res.append(complex(re[i], im[i]))
    return np.array(res)


def Initialise_V():
    V_List = [np.zeros((2, 1))]
    for lam in HALF_LAMBDA:
        e, U = eigh(F(lam))
        square_root_M = np.diag([cmath.sqrt(x) for x in e])
        if lam == np.pi:
            l = [np.random.normal(scale = np.sqrt(2)) for _ in np.arange(0, q)]
            V = np.matmul(np.matmul(U, square_root_M), np.array(l))
        else:
            V = np.matmul(np.matmul(U, square_root_M), Randomise_Z())
        V_List.append(V)
    reversed_V_List = reversed(V_List[1:-1])
    for item in reversed_V_List:
        V_List.append(np.conjugate(item))
    return V_List


def Weighted_Sum(V_List, t):
    return np.dot(np.exp(1j * t * LAMBDA)[1:-1], V_List[1:])
    # res = V_List[1] * np.exp(1j * t * LAMBDA[1])
    # for k in np.arange(2, n):
    #     res += V_List[k] * np.exp(1j * t * LAMBDA[k])
    # return res


def Generate_One_Path(V_List):
    X = np.zeros(shape = (2, n + 1))
    for t in np.arange(1, n + 1, 1):
        tmp = np.sqrt(np.pi / n) * Weighted_Sum(V_List, t)
        X[:, t] = tmp
        print(tmp)
    return X


# def argand(a):
#     import matplotlib.pyplot as plt
#     import numpy as np
#     for x in range(len(a)):
#         plt.plot([0, a[x].real], [0, a[x].imag], 'ro-', label = 'python')
#     limit = np.max(np.ceil(np.absolute(a)))  # set limits for axis
#     plt.xlim((-limit, limit))
#     plt.ylim((-limit, limit))
#     plt.ylabel('Imaginary')
#     plt.xlabel('Real')
#     plt.show()


if __name__ == '__main__':
    # print(Randomise_Z())
    V_List = Initialise_V()
    # for item in enumerate(V_List):
    #     print(item)
    X = Generate_One_Path(V_List)
    # print(X)

    # a = np.exp(-1j * (-np.pi + LAMBDA[1:]))
    # argand(a)
