#################################
# Your name: Jonathan Makovsky
#################################

import numpy as np
from matplotlib import pyplot as plt
from process_data import parse_data
from scipy.special import softmax

np.random.seed(7)


def adaboost(X_train, y_train, T):
    hypo_arr = []
    alphas_arr = []
    size = y_train.size
    theta_arr = np.unique(X_train)
    dt = distribution(size)

    for t in np.arange(T):
        h_pred, h_index, theta, ht = weak_learner(dt, X_train, y_train, theta_arr)
        hypo_arr.append((h_pred, h_index, theta))
        eps_t = err_prob(dt, y_train, ht)
        wt = (1/2)*np.log((1-eps_t)/eps_t)
        alphas_arr.append(wt)
        dt = softmax(np.log(dt)+(-1)*wt*y_train*ht)
    
    return (np.array(hypo_arr), np.array(alphas_arr))



def error_h(dt, X_train, y_train, theta, index):
    pos_err = 0
    neg_err = 0
    
    pos_err = err_prob(dt, y_train, ht(1, index, theta, X_train))
    neg_err = err_prob(dt, y_train, ht(-1, index, theta, X_train))

    return 2 * (pos_err<=neg_err)-1, min(pos_err, neg_err) 


def distribution(n):
    return np.full(n, 1/n)

def ht(pred, index, theta, X):
    return 2 * (X[:,int(index)]<=theta)-1 if pred==1 else 2*(X[:,int(index)]>theta)-1

def best_theta_for_index(dt, X_train, y_train, index, theta_arr):
    best_err = 1
    best_theta = -1
    best_pred = 0

    for theta in theta_arr:
        pred, err = error_h(dt, X_train, y_train, theta, index)
        if err < best_err:
            best_pred, best_theta, best_err = (pred, theta, err)
    return best_pred, best_theta, best_err

def final_error(X, y, hypo, alpha_arr, T, loss='zo'):
    classify_X_arr = classify_x(X, hypo, alpha_arr, T)
    if loss == 'zo':
        return np.array([np.average(classify_X_arr[t]!=y) for t in range(T)])
    elif loss == 'exp':
        mat = classify_X_arr * np.tile(alpha_arr, (classify_X_arr.shape[1], 1)).transpose()
        return np.array([np.average(np.exp(-1*y*np.sum(mat[:t+1,:], axis=0))) for t in range(T)])
    return


def err_prob(dt, y, ht):
    return np.sum(dt*(y!=ht))

def classify_x(X, hypo_arr, alpha_arr, T):
    summ = np.zeros(X.shape[0])
    classify_t_arr = []

    for t in np.arange(T):
        h_pred, h_index, theta = hypo_arr[t]
        summ += alpha_arr[t]*ht(h_pred, h_index, theta, X)
        classify_t_arr.append(2*(summ>=0)-1)
    
    return np.array(classify_t_arr)

def weak_learner(dt, X_train, y_train, theta_arr):
    h_pred, h_index, theta, best_err = (0, -1, -1, 1)

    for index in np.arange(y_train.size):
        pred, theta, err = best_theta_for_index(dt, X_train, y_train, index, theta_arr)
        if err < best_err:
            h_pred, h_index, theta, best_err = (pred, index, theta, err)

    ht = ht(h_pred, h_index, theta, X_train)
    
    return h_pred, h_index, theta, ht


    
def graph(title, x, y1, y2, label_x, label_y1, label_y2):
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel("Error")
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.legend({label_y1, label_y2})
    plt.show()

#Section A
def sectionA(X_train, y_train, X_test, y_test, hypotheses, alpha_arr, T):
    train_errors = final_error(X_train, y_train, hypotheses, alpha_arr, T)
    test_errors = final_error(X_test, y_test, hypotheses, alpha_arr, T)
    graph("Train and Test Error wrt t", np.arange(T), train_errors, test_errors,"T", "train error", "test error")

#Srction B 
def sectionB(hypotheses, T, vocab):
    for t in range(T):  
        h_pred, h_index, theta = hypotheses[t]
        print("Weak classifier t={0}:\th_pred={1}\th_index={2}[{3}]\ttheta={4}".format(t,h_pred,h_index,vocab[h_index],theta))

#Section C
def sectionC(X_train, y_train, X_test, y_test, hypotheses, alpha_arr, T):
    train_errors = final_error(X_train, y_train, hypotheses, alpha_arr, T, 'exp')
    test_errors = final_error(X_test, y_test, hypotheses, alpha_arr, T, 'exp')
    graph("Train and Test loss exp Error wrt t", np.arange(T), train_errors, test_errors,"T", "train error", "test error")




def main():
    T = 80
    data = parse_data()
    if not data:
        return
    (X_train, y_train, X_test, y_test, vocab) = data

    hypo, alpha_arr = adaboost(X_train, y_train, T)

    
    sectionA(X_train, y_train, X_test, y_test, hypo, alpha_arr, T)
    sectionB(hypo, 10, vocab)
    sectionC(X_train, y_train, X_test, y_test, hypo, alpha_arr, T)
    

if __name__ == '__main__':
    main()

