import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import sys
from matplotlib.pyplot import cm
import os, glob
import shutil

dirs = glob.glob('./output/Q1ManyObject/*')
dir_nums = np.array([int(i[i.rfind('/')+1:]) for i in dirs])
mask = dir_nums % 10 == 0

for i in range(len(dirs)):
    if mask[i] == 0:
        shutil.rmtree(dirs[i])



# Q1ManyMany = """20 12070 9.20671
# 40 13283 10.132
# 80 21738 16.5812 
# 160 11918 9.09077"""

# Q1FullScale = """400 2 0 0
# 1600 2 1 0.00460829
# 6400 2 2 0.0235294 
# 25600 7 40 0.162602  
# 102400 32 689 1.22163
# 409600 122 11276 8.60107"""


# Q2ManyMany = """20 31844 9.72
# 40 32842 10.02
# 80 29948 9.14
# 160 33304 10.17"""

# Q2FullScale = """3200 2 1 0.0128205
# 12800 7 11 0.0714286
# 51200 25 111 0.390845
# 204800 118 1729 2.39806
# 819200 469 5576 21.2015"""


# Q3ManyMany = """40 8089 26.96
# 80 7538 24.96
# 160 7188 23.88
# 320 7474 24.76
# 640 6723 22.26"""

# Q3FullScale = """3200 3 0 0
# 12800 18 2 0.0204163
# 51200 72 26 0.142041
# 204800 312 254 0.789555
# 819200 1296 2880 9.50495"""


# def plot_regressions(num_list, times_list, deg, color):

#     if deg == 1:
#         A = np.vstack([num_list, np.ones(len(num_list))]).T
#         coefs, res, rank, s = np.linalg.lstsq(A, times_list, rcond=None)
#         # a, b, c = coefs
#         c, d = coefs
#         x = np.linspace(min(num_list), max(num_list), 100)
#         plt.plot(x, c*x + d, color=color, label='Linear SSE = {0:.2e}, LC = {1:.2e})'.format(float(res), float(c)), linestyle='solid')
#     elif deg == 2:
#         A = np.vstack([num_list**2, num_list, np.ones(len(num_list))]).T
#         coefs, res, rank, s = np.linalg.lstsq(A, times_list, rcond=None)
#         # a, b, c = coefs
#         b, c, d = coefs
#         x = np.linspace(min(num_list), max(num_list), 100)
#         plt.plot(x, b*(x**2) + c*x + d, color=color, label='Quadratic (SSE = {0:.2e}, LC = {1:.2e})'.format(float(res), float(b)), linestyle='dotted')
#     elif deg == 3:
#         A = np.vstack([num_list**3, num_list**2, num_list, np.ones(len(num_list))]).T
#         coefs, res, rank, s = np.linalg.lstsq(A, times_list, rcond=None)
#         print(coefs)
#         a, b, c, d = coefs
#         print(a)
#         print(b)
#         print(c)
#         print(d)
#         print(res)
#         x = np.linspace(min(num_list), max(num_list), 100)
#         plt.plot(x, a*(x**3) + b*(x**2) + c*x + d, color=color, label='Cubic $SSE = 0, LC = {0:.2e}$)'.format(float(a)), linestyle='dashed')


# table = sys.argv[1]
# data = eval(table).split('\n')
# print(data)

# matrix = np.array([[float(j) for j in data[i].split()] for i in range(len(data))])
# print(matrix)


# def plot_basic_stuff():

#     color = cm.rainbow(np.linspace(0,1,10))
#     if matrix.shape[1] == 3:
#         # Scatter plot the data points
#         fig, ax = plt.subplots()
#         ax.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
#         ax.scatter(matrix[:,0], matrix[:,1])
#         plt.xlabel('Number of Objects')
#         plt.ylabel('Total time')
#         plt.title('Total Time vs Grid Size in Object Number Scaling Experiment')

#         for deg in range(1, 4):
#             plot_regressions(matrix[:,0], matrix[:,1], deg, color[deg-1])
        
#         plt.legend()
#         plt.savefig('TotalTime{0}Graph.png'.format(table))

#         plt.clf()

        
#         fig, ax = plt.subplots()
#         ax.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
#         ax.scatter(matrix[:,0], matrix[:,2])
#         plt.xlabel('Number of Objects')
#         plt.ylabel('Average time')
#         plt.title('Average Time vs Grid Size in Object Number Scaling Experiment')


#         for deg in range(1, 4):
#             plot_regressions(matrix[:,0], matrix[:,2], deg, color[deg-1])
        
#         plt.legend()

#         plt.savefig('AverageTime{0}Graph.png'.format(table))
#     else:

#         # Scatter plot the data points
#         fig, ax = plt.subplots()
#         ax.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
#         ax.scatter(matrix[:,0], matrix[:,2])
#         plt.xlabel('$N$')
#         plt.ylabel('CPU Time (s)')
#         plt.title('Total Time vs Grid Size in Full Scaling Experiment')

#         for deg in range(1, 4):
#             plot_regressions(matrix[:,0], matrix[:,2], deg, color[deg-1])
        
#         plt.legend()
#         plt.savefig('TotalTime{0}Graph.png'.format(table))

#         plt.clf()

        
#         fig, ax = plt.subplots()
#         ax.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
#         ax.scatter(matrix[:,0], matrix[:,3])
#         plt.xlabel('$N$')
#         plt.ylabel('Average CPU Time Per Step')
#         plt.title('Average Time vs Grid Size in Full Scaling Experiment')

#         for deg in range(1, 4):
#             plot_regressions(matrix[:,0], matrix[:,3], deg, color[deg-1])
        
#         plt.legend()

#         plt.savefig('AverageTime{0}Graph.png'.format(table))
        


#         # plt.xscale('log')
#         plt.show()


# # def plot_time_stuff():
# plot_basic_stuff()
