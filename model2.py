import numpy as np
from scipy.io import loadmat, savemat
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from joblib import dump, load
import h5py
import time

def model_predict():
	"""
	Load trained model and run prediction on saved data.

	Returns:
		y_pred: predicted values 
		t: time to run the prediction
	"""

	# load mat file with filename for model to load
	data = loadmat('model_data.mat')
	# load saved model from filename
	model = load(data['modelname'][0])
	# load saved data to run prediction on
	with h5py.File('test_data.mat', 'r') as data:
		X = np.array(data['X'])

	# start timing
	tic = time.perf_counter()
	X = np.transpose(X,(0,2,1))
	X = X[:,1:,:]

	# scale data
	xnorm = np.amax(np.linalg.norm(X,np.inf,axis=1),axis=1)
	X = X/xnorm[:,None,None]
	X = np.reshape(X,(np.shape(X)[0],-1))

	if len(model) == 3: 
		ralphie = model[0]
		ymin = model[1]
		ymax = model[2]
		y_pred = ralphie.predict(X)*(ymax-ymin)-ymin	# output is scaled
	else: 
		y_pred = model.predict(X)
	
	# finish timing
	toc = time.perf_counter()
	t = toc-tic

	return y_pred, t


def main():
	# run model prediction and save outputs to mat file
	y_pred, t = model_predict()
	savemat('predicted_data.mat',{'y_pred':y_pred,'time':t})

if __name__=='__main__':
    main()
