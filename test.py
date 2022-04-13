import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np
import os

def install():
    subprocess.check_call([sys.executable, "-m", "pip", "install","-v","--editable","."])

def read(file_name):
    return np.loadtxt(file_name)

def plot_data(f1, f2, c1 ,c2):
    y = read(f1)
    x = np.asarray([k for k in range(len(y))])
    y2 = read(f2)
    x2 = np.asarray([k for k in range(len(y2))])
    ff, (ax1,ax2) = plt.subplots(2,1)
    ax1.plot(x,y, c1)
    ax1.plot(x2,y2, c2)
    m = len(y)
    if len(y) > len(y2):
        m= len(y2)
    diff = []
    x3 = []
    for k in range(m):
        v = 0.0
        if y[k] != 0:
            v = 100*(y[k]-y2[k])/y[k]
        diff.append(v)
        x3.append(k)
    ax2.plot(x3, diff,c2)

try:   
    os.remove("quasielastic_test.tx")
    os.remove("quasielastic_test.python.lpt")
    os.remove("quasielastic_test2.t")
    os.remove("quasielastic_test.python2.lpt")
except:
    pass
install()
subprocess.call(["python", "quasielasticbayes\\test\\qldata_test.py"])



plot_data("quasielastic_test.tx", "quasielastic_test.python.lpt", "r", "g--")
#plot_data("quasielastic_test2.t", "quasielastic_test.python2.lpt", "r", "g--")
#plot_data("quasielastic_test.python.lpt", "quasielastic_test.python2.lpt", "r", "g--")


plt.show()






