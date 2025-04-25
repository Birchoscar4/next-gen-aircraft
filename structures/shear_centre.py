import numpy as np

sy = 1
f_height = 1
r_height = 1
c_height = 1
f_halfchord = 0.3
r_halfchord = 0.3
f_theta = np.arctan((c_height-f_height)/(2*f_halfchord))
r_theta = np.arctan((c_height-r_height)/(2*r_halfchord))
Ixx = 0.00001
A_cap = 0.01

t_spar = 0.002
t_skin = 0.002

# lengths of sections
L12 = f_height
L23 = np.sqrt(f_halfchord**2+(0.5*(c_height-f_height))**2)
L34 = np.sqrt(r_halfchord**2+(0.5*(c_height-r_height))**2)
L45 = r_height
L56 = L34
L61 = L23

# Areas of sections
A1 = (t_skin*L61)/2 + (t_spar*L12)/2 + A_cap
A2 = (t_spar*L12)/2 + (t_skin*L23)/2 + A_cap
A3 = (t_skin*L23)/2 + (t_skin*L34)/2
A4 = (t_skin*L34)/2 + (t_spar*L45)/2 + A_cap
A5 = (t_spar*L45)/2 + (t_skin*L56)/2 + A_cap
A6 = (t_skin*L56)/2 + (t_skin*L61)/2

A_tot = A1+A2+A3+A4+A5+A6

#y-pos of points
y1 = -f_height/2
y2 = f_height/2
y3 = c_height/2
y4 = r_height/2
y5 = -r_height/2
y6 = -c_height/2
Ycg=(A1*y1+A2*y2+A3*y3+A4*y4+A5*y5+A6*y6)/A_tot

#x-pos of points
x1 = -f_halfchord
x2 = -f_halfchord
x3 = 0
x4 = r_halfchord
x5 = r_halfchord
x6 = 0
Xcg=(A1*x1+A2*x2+A3*x3+A4*x4+A5*x5+A6*x6)/A_tot

# open shear flows
q1openx= 0
q2openx=q1openx-(sy/Ixx)*A2*(y2-Ycg)
q3openx=q2openx-(sy/Ixx)*A3*(y3-Ycg)
q4openx=q3openx-(sy/Ixx)*A4*(y4-Ycg)
q5openx=q4openx-(sy/Ixx)*A5*(y5-Ycg)
q6openx=q5openx-(sy/Ixx)*A6*(y6-Ycg)

# calculate p
a1 = - y2 + y1
b1 = x2 - x1
c1 = x1*-a1 - y1*b1
p1 = abs((a1*0+b1*0+c1)/np.sqrt(a1*a1 + b1*b1))
a2 = - y3 + y2
b2 = x3 - x2
c2 = x2*-a2 - y2*b2
p2 = abs((a2*0+b2*0+c2)/np.sqrt(a2*a2 + b2*b2))
a3 = - y4 + y3
b3 = x4 - x3
c3 = x3*-a3 - y3*b3
p3 = abs((a3*0+b3*0+c3)/np.sqrt(a3*a3 + b3*b3))
a4 = - y5 + y4
b4 = x5 - x4
c4 = x4*-a4 - y4*b4
p4 = abs((a4*0+b4*0+c4)/np.sqrt(a4*a4 + b4*b4))
a5 = - y6 + y5
b5 = x6 - x5
c5 = x5*-a5 - y5*b5
p5 = abs((a5*0+b5*0+c5)/np.sqrt(a5*a5 + b5*b5))
a6 = - y1 + y6
b6 = x1 - x6
c6 = x6*-a6 - y6*b6
p6 = abs((a6*0+b6*0+c6)/np.sqrt(a6*a6 + b6*b6))

enclosed_area = (f_height+c_height)*f_halfchord*0.5 + (r_height+c_height)*r_halfchord*0.5
print(enclosed_area)

q_s0 = -(q1openx*L12*p1+q2openx*L23*p2+q3openx*L34*p3+q4openx*L45*p4+q5openx*L56*p5+q6openx*L61*p6)/(2*enclosed_area)

q1 = q1openx + q_s0
q2 = q2openx + q_s0
q3 = q3openx + q_s0
q4 = q4openx + q_s0
q5 = q5openx + q_s0
q6 = q6openx + q_s0

# take moments about 0,0
x_shear = (q1*p1*L12+q2*p2*L23+q3*p3*L23+q4*p4*L45+q5*p5*L56+q6*p6*L61)/sy
print(x_shear)
