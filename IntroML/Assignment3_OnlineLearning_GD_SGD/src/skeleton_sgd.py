#################################
# Your name: Jonathan Makovsky
#################################


import numpy as np
import numpy.random
from sklearn.datasets import fetch_openml
import sklearn.preprocessing
import matplotlib.pyplot as plt
import scipy as sp
import math

"""
Please use the provided function signature for the SGD implementation.
Feel free to add functions and other code, and submit this file with the name sgd.py
"""

#------------------------------- Helper function difined by the course -------------------------------

def helper():
    mnist = fetch_openml('mnist_784', as_frame=False)
    data = mnist['data']
    labels = mnist['target']

    neg, pos = "0", "8"
    train_idx = numpy.random.RandomState(0).permutation(np.where((labels[:60000] == neg) | (labels[:60000] == pos))[0])
    test_idx = numpy.random.RandomState(0).permutation(np.where((labels[60000:] == neg) | (labels[60000:] == pos))[0])

    train_data_unscaled = data[train_idx[:6000], :].astype(float)
    train_labels = (labels[train_idx[:6000]] == pos) * 2 - 1

    validation_data_unscaled = data[train_idx[6000:], :].astype(float)
    validation_labels = (labels[train_idx[6000:]] == pos) * 2 - 1

    test_data_unscaled = data[60000 + test_idx, :].astype(float)
    test_labels = (labels[60000 + test_idx] == pos) * 2 - 1

    # Preprocessing
    train_data = sklearn.preprocessing.scale(train_data_unscaled, axis=0, with_std=False)
    validation_data = sklearn.preprocessing.scale(validation_data_unscaled, axis=0, with_std=False)
    test_data = sklearn.preprocessing.scale(test_data_unscaled, axis=0, with_std=False)
    return train_data, train_labels, validation_data, validation_labels, test_data, test_labels

# --------------------------------- End of helper function -------------------------------------

def SGD_hinge(data, labels, C, eta_0, T):
    W = np.zeros(data.shape[1])
    for t in range(1, T+1):
        W = update_w(W, t, data, labels, C, eta_0)
    return W


def SGD_log(data, labels, eta_0, T, arr = np.array([])):
    W = np.zeros(data.shape[1])
    counter = 0
    for t in range (1, T+1):
        W = update_w_loss(W, t, data, labels, eta_0)
        if (np.size(arr) > 0):
            arr[counter] = np.linalg.norm(W)
            counter += 1
    return W

#################################
#********************** Section1 ************************************
# --------------------- Function that updates W according to section 1

def update_w(W, t, data, labels, C, eta_0):
    eta_t = eta_0 / t
    index_arr = np.arange(data.shape[0])
    i = np.random.choice(index_arr)
    if ((W.dot(data[i])) * labels[i] < 1):
        W = (1 - eta_t) * W + (eta_t * C * (labels[i] * (data[i]))) 
    else:
        W = (1 - eta_t) * W
    return W

#--------------------- Section 1_a ---------------------------

def section1_a(train_data, train_labels, validation_data, validation_labels):
    iter = 10
    C = 1
    T = 1000
    values = np.logspace(-1, 1, 50)
    acc_array = np.zeros(values.size)
    counter = 0
    for eta_0 in values:   
        for i in range(iter):
            W = SGD_hinge(train_data,train_labels,C, eta_0, T)
            err = err_assessment(W, validation_data, validation_labels)
            acc_array[counter] = err 
        counter += 1
    
    plt.plot(values, acc_array)
    plt.xscale('log')
    plt.xlabel('eta_0')
    plt.ylabel('Accuracy')
    plt.title("Accuracy as a function of eta_0")
    plt.show()

    return values[np.argmax(acc_array)]


def err_assessment(W, validation_data, validation_labels):
    train_data = np.transpose(validation_data)
    return np.mean(W.dot(train_data) * (validation_labels) >= 0)
# -------------------- End of Section1_a ----------------------------
# -------------------- Section1_b -----------------------------------

def section1_b(train_data, train_labels, validation_data, validation_labels, eta_0):
    iter = 10
    T = 1000
    c_values = np.logspace(-4, -2, 50)
    acc_array = np.zeros(c_values.size)
    counter = 0
    for C in c_values:   
        for i in range(iter):
            W = SGD_hinge(train_data,train_labels,C, eta_0, T)
            err = err_assessment(W, validation_data, validation_labels)
            acc_array[counter] = err 
        counter += 1
    
    plt.plot(c_values, acc_array)
    plt.xscale('log')
    plt.xlabel('C')
    plt.ylabel('Accuracy')
    plt.title("Accuracy as a function of C")
    plt.show()

    return c_values[np.argmax(acc_array)]

# ------------------------- End of Section1_b -----------------------------

# ------------------------- Section1_c ------------------------------------

def section1_c(C, train_data, train_labels, validation_data, validation_labels, eta_0):
    T = 20000
    W = SGD_hinge(train_data, train_labels, C, eta_0, T)
    plt.imshow(np.reshape(W, (28, 28)), interpolation='nearest')
    plt.show()
    return W

# ------------------------- End of Srction1_c -----------------------------
#################################

train_data, train_labels, validation_data, validation_labels, test_data, test_labels = helper()

eta_0 = section1_a(train_data, train_labels, validation_data, validation_labels)
C = section1_b(train_data, train_labels, validation_data, validation_labels, eta_0)
w_res = section1_c(C, train_data, train_labels, validation_data, validation_labels, eta_0)

# Srction1_d
acc = err_assessment(w_res, test_data, test_labels)
print("The accuracy of the best classifier on the test set is: ", acc)


#********************************* End of Section1 ************************

#********************************* Section2 *******************************
def update_w_loss(W, t, data, labels, eta_0):
    eta_0 = eta_0 / t
    index_arr = np.arange(data.shape[0])
    i = np.random.choice(index_arr)
    exp = -labels[i] * (W.dot(data[i]))
    gradient = (-1*labels[i] * data[i]) * (sp.special.softmax([0, exp])[1]) 
    res = W - (gradient *(eta_0))
    return res
#--------------------------------- Section2_a -----------------------------
def section2_a(train_data, train_labels, validation_data, validation_labels):
    iter = 10
    T = 1000
    eta_0_array = np.logspace(-3, 1, 50)
    acc_array = np.zeros(eta_0_array.size)
    counter = 0
    for eta_0 in eta_0_array:   
        for i in range(iter):
            W = SGD_log(train_data,train_labels, eta_0, T)
            err = err_assessment(W, validation_data, validation_labels)
            acc_array[counter] = err 
        counter += 1
    
    plt.plot(eta_0_array, acc_array)
    plt.xscale('log')
    plt.xlabel('eta_0')
    plt.ylabel('Accuracy')
    plt.title("Accuracy as a function of eta_0")
    plt.show()

    return eta_0_array[np.argmax(acc_array)]

    
#--------------------------------- End of Section2_a ----------------------

#--------------------------------- Section2_b -----------------------------
def section_2b(train_data, train_labels, eta_0):
    T = 20000
    W = SGD_log(train_data, train_labels, eta_0, T)
    plt.imshow(np.reshape(W, (28, 28)), interpolation='nearest')
    plt.show()
    return W
#--------------------------------- End of Section2_b ----------------------

#--------------------------------- Section2_c -----------------------------
def section_2c(train_data, train_labels, eta_0, arr):
    T = 20000
    W = SGD_log(train_data, train_labels, eta_0, T, arr)
    return W

#--------------------------------- End of Section2_c ----------------------

#********************************* End of Section2 ************************
train_data, train_labels, validation_data, validation_labels, test_data, test_labels = helper()

eta_0 = section2_a(train_data, train_labels, validation_data, validation_labels)
#W = section_2b(train_data, train_labels, eta_0)
#acc = err_assessment(W, test_data, test_labels)
#print("The accuracy of the best classifier on the test set is: ", acc)
arr = np.zeros(20000)
W = section_2c(train_data, train_labels, eta_0, arr)
X = np.arange(20000)
plt.plot(X, arr)
#plt.xscale('log')
plt.xlabel('Iteration')
plt.ylabel('W norm')
plt.title("W norm as function of the iteration")
plt.show()

