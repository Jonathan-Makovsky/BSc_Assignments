import numpy as np
import numpy.random
from matplotlib import pyplot as plt
from sklearn.datasets import fetch_openml
from scipy.spatial.distance import euclidean

#Copied from the assignment
mnist = fetch_openml('mnist_784', as_frame=False)
data = mnist['data']
labels = mnist['target']
idx = numpy.random.RandomState(0).choice(70000, 11000)
train = data[idx[:10000], :].astype(int)
train_labels = labels[idx[:10000]]
test = data[idx[10000:], :].astype(int)
test_labels = labels[idx[10000:]]
#************

#K-NN algo

def label_prediction(images, labels, query, k):
    distances = []
    for curr in images:  
        distance = euclidean(curr, query)
        distances.append(distance)
    sort = np.argsort(distances)[:k]
    close_label = [labels[i] for i in sort]
    res = max(set(close_label), key = close_label.count) 
    return res


#Q2_b

def accuracy_test(test_images, test_labels, train_images, train_labels, k):
    c = 0
    size = len(test_images)
    for i in range(size):
        predict = label_prediction(train_images, train_labels, test_images[i], k)
        if (test_labels[i] == predict): 
            c += 1 
    res = c / size
    return res

acc = accuracy_test(test, test_labels, train[:1000], train_labels[:1000], k = 10)
print("The accuracy: ", acc)

#Q2_c
Y = np.zeros(100) 
x = np.linspace(1,100, 100) 
for i in range(1,101):
    Y[i - 1] = accuracy_test(test, test_labels, train[:1000], train_labels[:1000], k = i)

plt.title("prediction accuracu (by k)")
plt.xlabel("k")
plt.ylabel("Accuracy")
plt.plot(x,Y,'o')
plt.plot(x,Y)
plt.show()

#Q2_d
Y = np.zeros(50) 
x = np.linspace(100,5000, 50) 
for i in range(100,5001,100):
    Y[(i//100) - 1] = accuracy_test(test, test_labels,train[:i], train_labels[:i], k = 1)

plt.title("prediction accuracy (by n)")
plt.xlabel("n")
plt.ylabel("Accuracy") 
plt.plot(x,Y,'o')
plt.plot(x,Y)
plt.show()
#-----------------------------------------------------------------