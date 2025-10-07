import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import sympy as sy

def bernstein_poly(n, i, t):
    # 计算bernstein基多项式
    return comb(n, i) * (t ** i) * ((1 - t) ** (n - i))

def comb(n, k):
    # 计算组合数 C(n, k)
    from math import factorial
    return factorial(n) // (factorial(k) * factorial(n - k))

class Bezier3D:
    def __init__(self, control_points):
        self.control_points = np.array(control_points)
        self.n = len(control_points) - 1
        
    def point(self, t):
        # 计算Bezier曲线点
        p = np.zeros(3)
        for i in range(self.n + 1):
            b = bernstein_poly(self.n, i, t)
            p += b * self.control_points[i]
        return p
    
    def derivative(self, t):
        # 计算一阶导数
        p = np.zeros(3)
        n = self.n
        for i in range(n):
            b = bernstein_poly(n-1, i, t)
            p += b * n * (self.control_points[i+1] - self.control_points[i])
        return p
    
    def second_derivative(self, t):
        # 计算二阶导数
        p = np.zeros(3)
        n = self.n
        for i in range(n-1):
            b = bernstein_poly(n-2, i, t)
            p += b * n * (n-1) * (self.control_points[i+2] - 2*self.control_points[i+1] + self.control_points[i])
        return p

def curvature(der1, der2):
    # 计算曲率
    cross = np.linalg.norm(np.cross(der1, der2))
    denom = np.linalg.norm(der1) ** 3
    if denom == 0:
        return 0
    return cross / denom

def Show_Curvature(control_points):
    bezier = Bezier3D(control_points)
    ts = np.linspace(0, 1, 100)
    curvatures = []
    for t in ts:
        d1 = bezier.derivative(t)
        d2 = bezier.second_derivative(t)
        k = curvature(d1, d2)
        curvatures.append(k)

    # 绘制曲率曲线
    plt.plot(ts, curvatures)
    plt.xlabel('Parameter t')
    plt.ylabel('Curvature κ(t)')
    plt.title('3D Bezier Curve Curvature')
    plt.grid(True)
    plt.show()

def rotate_vector_around_axis(v, k, theta):
    k = k / np.linalg.norm(k)  # 单位化旋转轴
    v = np.array(v, dtype=float)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    v_rot = v * cos_theta + np.cross(k, v) * sin_theta + k * np.dot(k, v) * (1 - cos_theta)
    return v_rot

def GeneratePlaneClassA(N,vs,ve,T,ps,pe):
    vs = vs/np.linalg.norm(vs)
    ve = ve/np.linalg.norm(ve)
    T = T/np.linalg.norm(T)
    theta = np.arccos(np.dot(vs,ve))/(N-1)
    if(np.dot(np.cross(vs,ve),T)<0):
        theta=-theta
    X = (pe-ps)/np.linalg.norm(pe-ps)
    Y = np.cross(T,X)
    Y = Y/np.linalg.norm(Y)
    vectors = [vs]
    for i in range(1,N):
        vectors.append(rotate_vector_around_axis(vectors[-1],T,theta))
    s = sy.symbols('x')
    x_s = (np.dot(Y,vectors[-1]))*(s**(N-1))
    for i in range(N-1):
        x_s+=(np.dot(Y,vectors[i]))*(s**i)
    x_s_der = sy.diff(x_s)
    r = 1
    while(abs(x_s.subs(s,r))>1e-6):
        r = r-(x_s.subs(s,r))/(x_s_der.subs(s,r))
    print(r,abs(x_s.subs(s,r)))
    sum=0
    for i in range(N):
        sum+=np.dot(X,vectors[i])*(r**i)
    L = np.dot(pe-ps,X)/sum
    points = [ps]
    for i in range(N):
        points.append(points[-1]+vectors[i]*L*(r**i))
    points[-1]=pe
    return points,r

def GenerateSpaceClassA(NN,vs,ve,ps,pe):
    X = (vs+ve)/np.linalg.norm(vs+ve)
    Y = (np.cross(vs,ve))
    Y = Y/np.linalg.norm(Y)
    N = 30
    t_start = 0
    t_end = np.pi
    data_collect = []
    rep = []
    for iter in range(6):
        print('iter: ',iter)
        for i in range(1,N):
            t_theta = t_start*(1-i/N)+t_end*(i/N)
            T = np.cos(t_theta)*X+np.sin(t_theta)*Y
            vs_pro = vs-np.dot(vs,T)*T
            ve_pro = ve-np.dot(ve,T)*T
            vs_pro = vs_pro/np.linalg.norm(vs_pro)
            ve_pro = ve_pro/np.linalg.norm(ve_pro)
            H = np.dot(pe-ps,T)
            [cp,s] = GeneratePlaneClassA(NN,vs_pro,ve_pro,T,ps,pe-H*T)
            cp = [p.astype(np.float64) for p in cp]
            h = H*(s-1)/(s**NN-1)
            points=[ps]
            for k in range(1,NN+1):
                points.append(cp[k]+T*h*(s**k-1)/(s-1))
            points = [p.astype(np.float64) for p in points]
            re_vs = (points[1]-points[0])/np.linalg.norm(points[1]-points[0])
            re_ve = (points[-1]-points[-2])/np.linalg.norm(points[-1]-points[-2])
            error = np.linalg.norm(re_vs-vs)+np.linalg.norm(re_ve-ve)
            data_collect.append((t_theta,error,points))
        data_collect = sorted(data_collect,key=lambda x:x[1])
        print('t_theta: ',data_collect[0][0],' min_error: ',data_collect[0][1])
        t_start = data_collect[0][0]
        t_end = data_collect[1][0]
        rep = data_collect[0][2]
        data_collect.clear()
    return rep




vs = np.array([1,0,0])
ve = np.array([2,3,1])
ve = ve/np.linalg.norm(ve)
ps = np.array([0,0,0])
pe = np.array([10,2,1])
Controlp = GenerateSpaceClassA(6,vs,ve,ps,pe)
Show_Curvature(Controlp)