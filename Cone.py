import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import sympy as sy

epsilon = 0.000001

def NearlyEq(a,b):
    return (abs(a-b)<epsilon)

def plot_cone(axis_vector, height, angle_rad, ax=None):
    """
    绘制圆锥

    参数:
    axis_vector: np.array, 3维向量,圆锥的轴线方向
    height: float, 圆锥高度
    angle_rad: float, 轴线和母线夹角,单位弧度
    ax: matplotlib 3d坐标轴对象,可选
    """
    axis_vector = np.array(axis_vector)
    axis_vector = axis_vector / np.linalg.norm(axis_vector)  # 归一化轴向量

    # 圆锥顶点
    apex = np.array([0, 0, 0])

    # 底面圆心（沿轴线方向移动height距离）
    base_center = apex + axis_vector * height

    # 计算底面半径： r = height * tan(angle)
    radius = height * np.tan(angle_rad)

    # 构造垂直轴线的两个正交单位向量
    # 先找一个与轴线不共线的向量
    if (axis_vector == np.array([0,0,1])).all():
        not_colinear = np.array([1,0,0])
    else:
        not_colinear = np.array([0,0,1])

    v = np.cross(axis_vector, not_colinear)
    v = v / np.linalg.norm(v)
    w = np.cross(axis_vector, v)

    # 在底面坐标系中取一系列theta角度，画一个圆
    thetas = np.linspace(0, 2*np.pi, 100)
    circle_points = base_center[:, np.newaxis] + radius * (np.outer(v, np.cos(thetas)) + np.outer(w, np.sin(thetas)))

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    # 绘制圆锥侧面（三角面片）
    for i in range(len(thetas)-1):
        tri_x = [apex[0], circle_points[0, i], circle_points[0, i+1]]
        tri_y = [apex[1], circle_points[1, i], circle_points[1, i+1]]
        tri_z = [apex[2], circle_points[2, i], circle_points[2, i+1]]
        ax.plot_trisurf(tri_x, tri_y, tri_z, color='c', alpha=0.5)

    # 绘制底面圆
    ax.plot(circle_points[0, :], circle_points[1, :], circle_points[2, :], 'b')

    # 设置坐标轴比例相同
    max_range = max(height, radius) * 1.2
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

def ComputeAngle(vectors):
    re = []
    for i in range(1,len(vectors)):
        re.append(np.dot(vectors[i-1],vectors[i])/np.linalg.norm(vectors[i-1])/np.linalg.norm(vectors[i]))
    return re

def ComputeLengthRatio(points):
    re = []
    for i in range(2,len(points)):
        re.append(np.linalg.norm(points[i]-points[i-1])/np.linalg.norm(points[i-1]-points[i-2]))
    return re

def ComputeRotatingAngle(vectors,T,alpha):
    re = []
    k = (np.sin(alpha))**2
    for i in range(1,len(vectors)):
        theta = np.dot(vectors[i-1],vectors[i])/np.linalg.norm(vectors[i-1])/np.linalg.norm(vectors[i])
        sign=1
        if(np.dot(T,np.cross(vectors[i-1],vectors[i]))<0):
            sign=-1
        re.append(np.arccos(1-(1-np.cos(theta))/k)*sign)
    return re

def ComputeLength(points):
    re = []
    for i in range(1,len(points)):
        re.append(np.linalg.norm(points[i]-points[i-1]))
    return re

def solve_trig_equation(A, B, C):
    R = math.sqrt(A**2 + B**2)
    if abs(C) > R:
        # 无解，因为 |cos(x-α)| ≤ 1，而 C/R > 1 或 < -1
        return []
    alpha = math.atan2(B, A)
    val = C / R
    
    solutions = []
    # cos θ = val 有两组解 θ = acos(val), 2π - acos(val)
    theta1 = math.acos(val)
    theta2 = -theta1
    
    # x - α = θ  =>  x = θ + α
    x1 = theta1 + alpha
    x2 = theta2 + alpha
    
    # 统一到[0, 2π)
    x1 = x1 % (2 * math.pi)
    x2 = x2 % (2 * math.pi)
    
    solutions.append(x1)
    if x2 != x1:
        solutions.append(x2)
    
    return solutions

def IterationPoints(points,T,alpha):
    n = len(points)
    vectors=[]
    for i in range(1,len(points)):
        v = points[i]-points[i-1]
        vectors.append(v)
    vs = vectors[0]
    ve = vectors[-1]
    vs = vs/np.linalg.norm(vs)
    ve = ve/np.linalg.norm(ve)
    PHI = ComputeRotatingAngle(vectors,T,alpha)
    LENGTH = ComputeLength(points)
    X = np.array([1,0,0])
    Y = np.cross(T,X)
    Y = Y/np.linalg.norm(Y)
    X = np.cross(Y,T)
    newpoints = []
    for p in points:
        newpoints.append(p)
    for i in range(2,n-2):
        v_2 = points[i+1]-points[i-1]
        L_v_2 = np.linalg.norm(v_2)
        v_2 = v_2/L_v_2
        phi = (PHI[i-2]+PHI[i-1]+PHI[i])/3
        s1 = LENGTH[i-1]/LENGTH[i-2]
        s2 = LENGTH[i]/LENGTH[i-1]
        s3 = LENGTH[i+1]/LENGTH[i]
        s = (s1+s2+s3)/3
        A = (1+s)**2-(np.dot(v_2,T)**2)*2*s*(1-np.cos(phi))
        B = (np.dot(v_2,T)**2)*(1+s**2+2*s*np.cos(phi))
        cos_a = np.sqrt(B/A)
        sin_a = np.sqrt(1-B/A)
        cos_theta = 1-(1-np.cos(phi))*(1-B/A)
        eq1_A = 1+s*np.cos(phi)
        eq1_B = -s*np.sin(phi)
        eq1_C = np.sqrt(1+s**2+2*s*cos_theta)/sin_a*np.dot(v_2,X)
        solve1 = (solve_trig_equation(eq1_A,eq1_B,eq1_C))
        eq2_A = s*np.sin(phi)
        eq2_B = 1+s*np.cos(phi)
        eq2_C = np.sqrt(1+s**2+2*s*cos_theta)/sin_a*np.dot(v_2,Y)
        solve2 = (solve_trig_equation(eq2_A,eq2_B,eq2_C))
        xi = 100
        for so1 in solve1:
            for so2 in solve2:
                if(NearlyEq(so1,so2)):
                    xi = so1
        if(xi==100):
            print("Oh!Fuck!!!!!!!!!!!")
        v_0 = cos_a* T+sin_a * (np.cos(xi)*X+np.sin(xi)*Y)
        #v_1 = cos_a* T+sin_a * (np.cos(xi+phi)*X+np.sin(xi+phi)*Y)
        VV = (L_v_2/np.sqrt(1+s**2+2*s*cos_theta))*v_0
        newpoints[i] = points[i-1]+VV
    all_length = ComputeLength(newpoints)
    sum = 0
    for i in range(1,len(all_length)):
        sum = sum + (all_length[i]/all_length[i-1])
    sum = sum/(len(all_length)-1)
    newpoints[1] = newpoints[0]+vs*all_length[1]/sum
    newpoints[-2] = newpoints[-1]-ve*all_length[-2]*sum
    return newpoints

def GenerateAxis(control_p):
    vectors=[]
    for i in range(1,len(control_p)):
        v = control_p[i]-control_p[i-1]
        v = v/np.linalg.norm(v)
        vectors.append(v)
    N = len(vectors)
    vs = vectors[0]
    ve = vectors[N-1]
    X = vs+ve
    X = X/np.linalg.norm(X)
    Y = np.cross(vs,ve)
    Y = Y/np.linalg.norm(Y)
    A=0
    B=0
    C=0
    for i in range(1,N-1):
        A = A + (np.dot(X,vectors[i]-vs))**2
        B = B + (np.dot(Y,vectors[i]-vs))**2
        C = C + 2 * np.dot(X,vectors[i]-vs) * np.dot(Y,vectors[i]-vs)
    target_min = (A+B)/2-np.sqrt(((A-B)/2)**2+(C/2)**2)
    phi = np.arctan((A-B)/C)
    if(C<0):
        phi = phi+np.pi
    target_theta = np.pi*0.75-phi*0.5
    T = X*np.cos(target_theta)+Y*np.sin(target_theta)
    cos_alpha = np.dot(vs,T)
    alpha = np.arccos(cos_alpha)
    print('min ' + str(target_min))
    if(alpha>np.pi/2):
        alpha = np.pi-alpha
        T = -T
    print("alpha ",alpha,"T ",T)
    return vectors,T,alpha

vs = np.array([1,0,0])
ve = np.array([2,3,1])
ve = ve/np.linalg.norm(ve)
ps = np.array([0,0,0])
pe = np.array([10,2,1])
X = (vs+ve)/np.linalg.norm(vs+ve)
Y = (np.cross(vs,ve))
Y = Y/np.linalg.norm(Y)
N = 30
for i in range(N - 1):
    t_theta = np.pi/N*(i+1)
    T = np.cos(t_theta)*X+np.sin(t_theta)*Y
    vs_pro = vs-np.dot(vs,T)*T
    ve_pro = ve-np.dot(ve,T)*T
    vs_pro = vs_pro/np.linalg.norm(vs_pro)
    ve_pro = ve_pro/np.linalg.norm(ve_pro)
    vm_pro = vs_pro+ve_pro
    vm_pro = vm_pro/np.linalg.norm(vm_pro)
    cos_alpha = np.dot(T,vs)
    vm = cos_alpha*T+np.sqrt(1-cos_alpha*cos_alpha)*vm_pro
    [a,b,c]=(pe-ps)/np.linalg.norm(pe-ps)
    s = sy.symbols('x')
    x_s = ve[0]*(s**2)+vm[0]*s+vs[0]
    y_s = ve[1]*(s**2)+vm[1]*s+vs[1]
    z_s = ve[2]*(s**2)+vm[2]*s+vs[2]
    target_func = (a*y_s-b*x_s)**2+(c*y_s-b*z_s)**2+(c*x_s-a*z_s)**2
    func_der = sy.diff(target_func)
    roots = sy.solve(func_der, s)
    real_roots = [sy.re(r) for r in roots if (NearlyEq(sy.im(r), 0) and sy.re(r)>0)]
    value = [target_func.subs(s, r).evalf() for r in real_roots]
    computed_v=np.array([([x_s.subs(s,r),y_s.subs(s,r),z_s.subs(s,r)]) for r in real_roots])
    computed_v = np.array(computed_v,dtype=float)
    for i in range(len(computed_v)):
        computed_v[i] = computed_v[i]/np.linalg.norm(computed_v[i])
    print("Root: ",real_roots," value: ",value," axis: ",T)
    print("Computed: ",computed_v," Real: ",[a,b,c])

