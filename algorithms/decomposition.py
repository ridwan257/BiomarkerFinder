import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler



class PLSClassifier(ClassifierMixin, BaseEstimator):
    def __init__(self, n_components=2, threshold=0.5, std=True, vip=True):
        super().__init__()
        self.n_components = n_components
        self.threshold = threshold
        self.std = std
        self.vip = vip
        self._pls = PLSRegression(n_components=self.n_components, scale=False)
        self._scaler = StandardScaler()
    
    def fit(self, X, y):
        self.classes_ = np.unique(y)
        if self.std:
            X_scaled = self._scaler.fit_transform(X)
            self._pls.fit(X_scaled, y)
        else:
            self._pls.fit(X, y)

        # exposing all attributes
        for attr in dir(self._pls):
            if not attr.startswith("_") and not callable(getattr(self._pls, attr)):  
                setattr(self, attr, getattr(self._pls, attr))
        
        if self.vip:
            ## Formula found from —
            ## Farrés, M., Platikanov, S., Tsakovski, S., and Tauler, R. (2015) 
            ## Comparison of the variable importance in projection (VIP) and of the selectivity ratio (SR) methods 
            ## for variable selection and interpretation. J. Chemometrics, 29: 528–536. doi: 10.1002/cem.2736.
            
            LVs = self._pls.x_scores_ 
            w = self._pls.x_weights_
            b_f = self._pls.y_loadings_.flatten()
            # # y_scaler = StandardScaler().fit_transform(y.reshape(len(y), 1)).ravel()
            # # SS = np.dot((LVs.T)**2, y)
            # # ssy = np.array([np.var(model.y_scores_[:, a]) for a in range(a)])
            SS = (b_f**2) *  np.sum(LVs**2, axis=0)
            
            SS_w_array = np.dot(w**2, SS)
            # # SS_w_array = np.sum(SS * w**2, axis=1)
            VIPi = X.shape[1] * (SS_w_array / np.sum(SS))
            self.vip_ = np.sqrt(VIPi)
        
        return self
    
    def predict(self, X):
        if self.std:
            X_scaled = self._scaler.transform(X)
            y_pred = self._pls.predict(X_scaled).ravel()
        else:
            y_pred = self._pls.predict(X).ravel()
         
        return np.where(y_pred > self.threshold, self.classes_[1], self.classes_[0])
    


    def predict_proba(self, X):
        y_pred = self._pls.predict(X).ravel()
        return np.column_stack((1 - y_pred, y_pred))