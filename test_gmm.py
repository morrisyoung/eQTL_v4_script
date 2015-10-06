import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn import mixture





if __name__ == "__main__":





	##============================= part#0: get the corr_list ==============================
	corr_list = []
	## fit a Gaussian into the data, "../result_init/para_init_train_cis_corr.txt"
	file = open("../result_init/para_init_train_cis_corr.txt", 'r')
	#file = open("../result_init/para_init_train_cis_corr_tissues/para_corr_1.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		corr = float(line[1])
		if np.isnan(corr):  # NOTE: there are some NaN
			continue
		corr_list.append(corr)
	file.close()
	corr_list = np.array(corr_list)





	##============================= part#1: PDF of the original data ==============================
	plt.hist(corr_list, histtype='step', bins=500, normed=1, alpha=0.5)
	plt.plot([], [], color='blue', label="original pdf")
	#plt.title("Hidden batch distribution after one sample (with the initialization we calculated)")
	#plt.xlabel("Value of hidden batch")
	#plt.ylabel("Frequency")
	#plt.axis([-1, 1, 0, 7])  # un-normalized version: y ~ [0, 300]; normalized version: y ~ [0, 7]
	#plt.show()








	##============================= part#2: get the two Gaussian mixture ==============================
	np.random.seed(1)
	g = mixture.GMM(n_components=2)

	##=========== transform the corr_list into the required format (a wired format) ===========
	list = []
	for corr in corr_list:
		list.append([corr])
	corr_list = np.array(list)

	obs = corr_list
	#obs = np.concatenate((np.random.randn(100, 1), 10 + np.random.randn(300, 1)))
	g.fit(obs)
	#GMM(covariance_type='diag', init_params='wmc', min_covar=0.001, n_components=2, n_init=1, n_iter=100, params='wmc', random_state=None, thresh=None, tol=0.001)


	##=========== get model parameters ===========
	## show the learned model:
	# weight:
	#array = np.round(g.weights_, 2)
	array = g.weights_
	#array([ 0.75,  0.25])
	weight1 = array[0]
	weight2 = array[1]
	print weight1,
	print weight2


	# mean:
	#array = np.round(g.means_, 2)
	array = g.means_
	#array([[ 10.05], [ 0.06]])
	mean1 = array[0][0]
	mean2 = array[1][0]
	print mean1,
	print mean2


	# var:
	##array = np.round(g.covars_, 2)
	array = g.covars_
	#array([[[ 1.02]], [[ 0.96]]])
	var1 = array[0][0]
	var2 = array[1][0]
	print var1,
	print var2



	##============================= part#3: plot the learned two Gaussian ==============================
	## plotting a Gaussian, determined by "mean", "variance" and "weight"
	mean = mean1
	variance = var1
	sigma = math.sqrt(variance)
	weight = weight1
	x = np.linspace(-1,1,100)
	plt.plot(x, weight * mlab.normpdf(x, mean, sigma), 'red')
	plt.plot([], [], color='red', label="this is one learned Gaussian")

	mean = mean2
	variance = var2
	sigma = math.sqrt(variance)
	weight = weight2
	x = np.linspace(-1,1,100)
	plt.plot(x, weight * mlab.normpdf(x, mean, sigma), 'green')
	plt.plot([], [], color='green', label="this is one learned Gaussian")



	plt.title("Distribution of Pearson correlation values, and the Gaussian mixture")
	plt.xlabel("Pearson correlation value")
	plt.ylabel("pdf")
	plt.axis([-1, 1, 0, 7])  # un-normalized version: y ~ [0, 300]; normalized version: y ~ [0, 7]
	plt.legend()
	plt.show()










	'''
	##============================= part#2: get the one Gaussian mixture ==============================
	np.random.seed(1)
	g = mixture.GMM(n_components=1)

	##=========== transform the corr_list into the required format (a wired format) ===========
	list = []
	for corr in corr_list:
		list.append([corr])
	corr_list = np.array(list)

	obs = corr_list
	#obs = np.concatenate((np.random.randn(100, 1), 10 + np.random.randn(300, 1)))
	g.fit(obs)
	#GMM(covariance_type='diag', init_params='wmc', min_covar=0.001, n_components=2, n_init=1, n_iter=100, params='wmc', random_state=None, thresh=None, tol=0.001)


	##=========== get model parameters ===========
	## show the learned model:
	# weight:
	#array = np.round(g.weights_, 2)
	array = g.weights_
	#array([ 0.75,  0.25])
	weight1 = array[0]
	print weight1

	# mean:
	#array = np.round(g.means_, 2)
	array = g.means_
	#array([[ 10.05], [ 0.06]])
	mean1 = array[0][0]
	print mean1

	# var:
	##array = np.round(g.covars_, 2)
	array = g.covars_
	#array([[[ 1.02]], [[ 0.96]]])
	var1 = array[0][0]
	print var1


	##============================= part#3: plot the learned two Gaussian ==============================
	## plotting a Gaussian, determined by "mean", "variance" and "weight"
	mean = mean1
	variance = var1
	sigma = math.sqrt(variance)
	weight = weight1
	x = np.linspace(-1,1,100)
	plt.plot(x, weight * mlab.normpdf(x, mean, sigma), 'red')
	plt.plot([], [], color='red', label="this is one learned Gaussian")



	plt.title("Distribution of Pearson correlation values, and the Gaussian mixture")
	plt.xlabel("Pearson correlation value")
	plt.ylabel("pdf")
	plt.axis([-1, 1, 0, 7])  # un-normalized version: y ~ [0, 300]; normalized version: y ~ [0, 7]
	plt.legend()
	plt.show()
	'''
