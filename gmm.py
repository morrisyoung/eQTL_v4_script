import numpy as np
from sklearn import mixture
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab






if __name__ == "__main__":




	corr_list = []
	## fit a Gaussian into the data, "../result_init/para_init_train_cis_corr.txt"
	file = open("../result_init/para_init_train_cis_corr.txt", 'r')
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




	'''
	##============================= part#1: PDF of the original data ==============================
	plt.hist(corr_list, histtype='step', bins=500)
	plt.title("Hidden batch distribution after one sample (with the initialization we calculated)")
	plt.xlabel("Value of hidden batch")
	plt.ylabel("Frequency")
	plt.axis([-1, 1, 0, 300])
	plt.show()
	'''




	##============================= part#2: Gaussian plotting ==============================
	## plotting a Gaussian, determined by "mean", "variance" and "weight"
	mean = 0
	variance = 0.1
	sigma = math.sqrt(variance)
	weight = 2
	x = np.linspace(-1,1,100)
	plt.plot(x, weight * mlab.normpdf(x, mean, sigma), 'r')
	plt.plot([], [], color='red', label="this is standard Gaussian")

	plt.legend()
	plt.show()







